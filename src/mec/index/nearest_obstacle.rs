//! Nearest-obstacle distance queries, accelerated by a bounding-volume
//! hierarchy.
//!
//! A uniform grid looks attractive and is wrong for this problem. The solver
//! asks for the distance at points deep between obstacles, where the nearest
//! obstacle is far away, so a grid walk visits a number of cells that grows with
//! the primitive count: measured on a regular polygon's centre in the sibling
//! crate, one query cost 57 ns at 16 segments and 48.8 us at 4096. The
//! hierarchy below descends only into nodes whose box can still hold a closer
//! primitive, so the cost tracks tree depth rather than table size.
//!
//! Points and segments share one hierarchy: a point becomes a degenerate box,
//! which keeps the traversal uniform and costs one comparison per point leaf.
use super::super::input::{ObstacleSet, SegmentIndex};

/// A handful of primitives per leaf keeps the tree shallow and lets a leaf test
/// finish within one cache line's worth of work.
const LEAF_SIZE: usize = 4;

/// Below this many primitives a linear scan beats a tree descent.
const LINEAR_SCAN_MAX: usize = 32;

#[derive(Debug, Clone, Copy)]
struct Node {
    min_x: f64,
    min_y: f64,
    max_x: f64,
    max_y: f64,
    /// Left child index for an internal node, or the first slot of the leaf run.
    left: u32,
    /// Right child index for an internal node; unused for a leaf.
    right: u32,
    /// Primitive count for a leaf. Zero marks an internal node.
    count: u32,
}

/// One leaf entry: either a point or a segment of the obstacle set.
#[derive(Debug, Clone, Copy)]
struct Prim(u32);

/// Nearest-obstacle distance queries over an obstacle set.
#[derive(Debug, Clone)]
pub struct NearestObstacleIndex {
    points: Vec<[f64; 2]>,
    segments: SegmentIndex,
    nodes: Vec<Node>,
    order: Vec<Prim>,
}

impl NearestObstacleIndex {
    /// Build the hierarchy over an obstacle set.
    pub fn new(set: &ObstacleSet) -> Self {
        let n_points = set.points.len();
        let n_segments = set.segments.len();
        let total = n_points + n_segments;

        if total == 0 {
            return Self {
                points: set.points.clone(),
                segments: set.segments.clone(),
                nodes: Vec::new(),
                order: Vec::new(),
            };
        }

        let boxes: Vec<(f64, f64, f64, f64)> = set
            .points
            .iter()
            .map(|p| (p[0], p[1], p[0], p[1]))
            .chain((0..n_segments).map(|i| {
                (
                    set.segments.bbox_minx[i],
                    set.segments.bbox_miny[i],
                    set.segments.bbox_maxx[i],
                    set.segments.bbox_maxy[i],
                )
            }))
            .collect();

        let mut order: Vec<Prim> = (0..total as u32).map(Prim).collect();

        let mut nodes = Vec::with_capacity(2 * total.div_ceil(LEAF_SIZE) + 1);
        build(&mut nodes, &mut order, &boxes, 0, total);

        Self {
            points: set.points.clone(),
            segments: set.segments.clone(),
            nodes,
            order,
        }
    }

    /// The segment table the index was built from.
    #[inline]
    pub fn segments(&self) -> &SegmentIndex {
        &self.segments
    }

    /// Squared distance from (x, y) to the closest obstacle and its ordinal.
    ///
    /// The ordinal numbers points first, then segments, matching the order used
    /// to build the index. Callers that only need the distance ignore it.
    #[inline]
    pub fn nearest_distance_sq(&self, x: f64, y: f64) -> Option<(f64, usize)> {
        let total = self.points.len() + self.segments.len();
        if total == 0 {
            return None;
        }
        if total <= LINEAR_SCAN_MAX {
            return Some(self.linear_scan(x, y));
        }
        let mut best_sq = f64::INFINITY;
        let mut best_idx = 0usize;
        self.descend(0, x, y, &mut best_sq, &mut best_idx);
        Some((best_sq, best_idx))
    }

    /// Nearest distance, or `None` once a primitive is proven within
    /// `threshold_sq`.
    ///
    /// A branch-and-bound caller often does not need the exact distance: it only
    /// needs to know whether a point can beat an incumbent. Once one primitive
    /// is found closer than that bar, the point is useless. `seed` is a nearby
    /// primitive ordinal to test first, typically the support of the parent cell.
    #[inline]
    pub fn nearest_distance_sq_above(
        &self,
        x: f64,
        y: f64,
        threshold_sq: f64,
        seed: Option<usize>,
    ) -> Option<(f64, usize)> {
        let total = self.points.len() + self.segments.len();
        if total == 0 {
            return None;
        }
        if total <= LINEAR_SCAN_MAX {
            let mut best_sq = f64::INFINITY;
            let mut best_idx = 0usize;
            for i in 0..total {
                let d_sq = self.prim_distance_sq(i, x, y);
                if d_sq <= threshold_sq {
                    return None;
                }
                if d_sq < best_sq {
                    best_sq = d_sq;
                    best_idx = i;
                }
            }
            return Some((best_sq, best_idx));
        }

        let mut best_sq = f64::INFINITY;
        let mut best_idx = 0usize;
        if let Some(si) = seed {
            if si < total {
                let d_sq = self.prim_distance_sq(si, x, y);
                if d_sq <= threshold_sq {
                    return None;
                }
                best_sq = d_sq;
                best_idx = si;
            }
        }
        if self.descend_above(0, x, y, &mut best_sq, &mut best_idx, threshold_sq) {
            return None;
        }
        Some((best_sq, best_idx))
    }

    #[inline]
    fn prim_distance_sq(&self, ordinal: usize, x: f64, y: f64) -> f64 {
        let n_points = self.points.len();
        if ordinal < n_points {
            let p = self.points[ordinal];
            let dx = x - p[0];
            let dy = y - p[1];
            dx * dx + dy * dy
        } else {
            self.segments
                .point_segment_distance_sq(ordinal - n_points, x, y)
        }
    }

    fn linear_scan(&self, x: f64, y: f64) -> (f64, usize) {
        let total = self.points.len() + self.segments.len();
        let mut best_sq = f64::INFINITY;
        let mut best_idx = 0usize;
        for i in 0..total {
            let d_sq = self.prim_distance_sq(i, x, y);
            if d_sq < best_sq {
                best_sq = d_sq;
                best_idx = i;
            }
        }
        (best_sq, best_idx)
    }

    fn descend(&self, node_idx: usize, x: f64, y: f64, best_sq: &mut f64, best_idx: &mut usize) {
        let node = self.nodes[node_idx];
        if bbox_dist_sq(&node, x, y) >= *best_sq {
            return;
        }
        if node.count > 0 {
            let start = node.left as usize;
            for prim in &self.order[start..start + node.count as usize] {
                let (ordinal, d_sq) = self.eval(prim, x, y);
                if d_sq < *best_sq {
                    *best_sq = d_sq;
                    *best_idx = ordinal;
                }
            }
            return;
        }
        // Visit the nearer child first: its result tightens the bound before the
        // far child is tested, which is what makes the pruning effective.
        let left = node.left as usize;
        let right = node.right as usize;
        let dl = bbox_dist_sq(&self.nodes[left], x, y);
        let dr = bbox_dist_sq(&self.nodes[right], x, y);
        let (first, second) = if dl <= dr {
            (left, right)
        } else {
            (right, left)
        };
        self.descend(first, x, y, best_sq, best_idx);
        if bbox_dist_sq(&self.nodes[second], x, y) < *best_sq {
            self.descend(second, x, y, best_sq, best_idx);
        }
    }

    fn descend_above(
        &self,
        node_idx: usize,
        x: f64,
        y: f64,
        best_sq: &mut f64,
        best_idx: &mut usize,
        threshold_sq: f64,
    ) -> bool {
        let node = self.nodes[node_idx];
        if bbox_dist_sq(&node, x, y) >= *best_sq {
            return false;
        }
        if node.count > 0 {
            let start = node.left as usize;
            for prim in &self.order[start..start + node.count as usize] {
                let (ordinal, d_sq) = self.eval(prim, x, y);
                if d_sq <= threshold_sq {
                    return true;
                }
                if d_sq < *best_sq {
                    *best_sq = d_sq;
                    *best_idx = ordinal;
                }
            }
            return false;
        }
        let left = node.left as usize;
        let right = node.right as usize;
        let dl = bbox_dist_sq(&self.nodes[left], x, y);
        let dr = bbox_dist_sq(&self.nodes[right], x, y);
        let (first, second) = if dl <= dr {
            (left, right)
        } else {
            (right, left)
        };
        if self.descend_above(first, x, y, best_sq, best_idx, threshold_sq) {
            return true;
        }
        if bbox_dist_sq(&self.nodes[second], x, y) < *best_sq
            && self.descend_above(second, x, y, best_sq, best_idx, threshold_sq)
        {
            return true;
        }
        false
    }

    #[inline]
    fn eval(&self, prim: &Prim, x: f64, y: f64) -> (usize, f64) {
        let ordinal = prim.0 as usize;
        (ordinal, self.prim_distance_sq(ordinal, x, y))
    }

    /// Ordinal and squared distance of every primitive within `max_dist_sq`.
    pub fn supporting_primitives(&self, x: f64, y: f64, max_dist_sq: f64) -> Vec<usize> {
        let mut out = Vec::new();
        let total = self.points.len() + self.segments.len();
        if total == 0 {
            return out;
        }
        if total <= LINEAR_SCAN_MAX {
            for i in 0..total {
                if self.prim_distance_sq(i, x, y) <= max_dist_sq {
                    out.push(i);
                }
            }
            return out;
        }
        self.collect(0, x, y, max_dist_sq, &mut out);
        out
    }

    fn collect(&self, node_idx: usize, x: f64, y: f64, max_dist_sq: f64, out: &mut Vec<usize>) {
        let node = self.nodes[node_idx];
        if bbox_dist_sq(&node, x, y) > max_dist_sq {
            return;
        }
        if node.count > 0 {
            let start = node.left as usize;
            for prim in &self.order[start..start + node.count as usize] {
                let (ordinal, d_sq) = self.eval(prim, x, y);
                if d_sq <= max_dist_sq {
                    out.push(ordinal);
                }
            }
            return;
        }
        self.collect(node.left as usize, x, y, max_dist_sq, out);
        self.collect(node.right as usize, x, y, max_dist_sq, out);
    }
}

/// Build the subtree owning `order[start..end]`, returning its node index.
///
/// Splitting on the widest axis of the primitive centres keeps the boxes tight,
/// which is what lets the bound prune.
fn build(
    nodes: &mut Vec<Node>,
    order: &mut [Prim],
    boxes: &[(f64, f64, f64, f64)],
    start: usize,
    end: usize,
) -> usize {
    let (min_x, min_y, max_x, max_y) = bounds_of(boxes, order, start, end);

    let idx = nodes.len();
    nodes.push(Node {
        min_x,
        min_y,
        max_x,
        max_y,
        left: 0,
        right: 0,
        count: 0,
    });

    let count = end - start;
    if count <= LEAF_SIZE {
        nodes[idx].left = start as u32;
        nodes[idx].count = count as u32;
        return idx;
    }

    let (mut cmin_x, mut cmin_y) = (f64::INFINITY, f64::INFINITY);
    let (mut cmax_x, mut cmax_y) = (f64::NEG_INFINITY, f64::NEG_INFINITY);
    for prim in &order[start..end] {
        let b = boxes[prim_index(prim)];
        let cx = (b.0 + b.2) * 0.5;
        let cy = (b.1 + b.3) * 0.5;
        cmin_x = cmin_x.min(cx);
        cmin_y = cmin_y.min(cy);
        cmax_x = cmax_x.max(cx);
        cmax_y = cmax_y.max(cy);
    }
    let axis = if cmax_x - cmin_x >= cmax_y - cmin_y {
        0
    } else {
        1
    };

    let mid = start + count / 2;
    order[start..end].select_nth_unstable_by(count / 2, |a, b| {
        let ca = centre(boxes[prim_index(a)], axis);
        let cb = centre(boxes[prim_index(b)], axis);
        ca.partial_cmp(&cb).unwrap_or(std::cmp::Ordering::Equal)
    });

    nodes[idx].left = nodes.len() as u32;
    nodes[idx].count = 0;
    build(nodes, order, boxes, start, mid);
    let right = build(nodes, order, boxes, mid, end);
    nodes[idx].right = right as u32;
    idx
}

#[inline]
fn prim_index(prim: &Prim) -> usize {
    prim.0 as usize
}

fn bounds_of(
    boxes: &[(f64, f64, f64, f64)],
    order: &[Prim],
    start: usize,
    end: usize,
) -> (f64, f64, f64, f64) {
    let mut min_x = f64::INFINITY;
    let mut min_y = f64::INFINITY;
    let mut max_x = f64::NEG_INFINITY;
    let mut max_y = f64::NEG_INFINITY;
    for prim in &order[start..end] {
        let b = boxes[prim_index(prim)];
        min_x = min_x.min(b.0);
        min_y = min_y.min(b.1);
        max_x = max_x.max(b.2);
        max_y = max_y.max(b.3);
    }
    (min_x, min_y, max_x, max_y)
}

#[inline]
fn centre(b: (f64, f64, f64, f64), axis: usize) -> f64 {
    if axis == 0 {
        (b.0 + b.2) * 0.5
    } else {
        (b.1 + b.3) * 0.5
    }
}

#[inline]
fn bbox_dist_sq(node: &Node, x: f64, y: f64) -> f64 {
    let dx = if x < node.min_x {
        node.min_x - x
    } else if x > node.max_x {
        x - node.max_x
    } else {
        0.0
    };
    let dy = if y < node.min_y {
        node.min_y - y
    } else if y > node.max_y {
        y - node.max_y
    } else {
        0.0
    };
    dx * dx + dy * dy
}
