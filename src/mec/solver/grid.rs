//! Grid solver: tiles the search box with cells and refines the best one.
//!
//! The search is a branch and bound on one invariant: the distance to the
//! nearest obstacle is 1-Lipschitz, so for a cell of half-size `h` centred at
//! `c`, no point of that cell can be farther from an obstacle than
//! `obstacle_distance(c) + h * sqrt(2)`. A cell leaves the queue only when its
//! bound proves it cannot beat the incumbent by more than the tolerance.
//!
//! Two design points separate the empty circle from the inscribed circle.
//!
//! **Region membership is a filter, not part of the bound.** A cell outside the
//! region can still hold an admissible point near its edge, so it stays in the
//! queue with a valid Lipschitz bound and prunes once that bound falls behind the
//! incumbent. Marking it hopeless up front would discard the boundary optima,
//! which is exactly where a constrained answer often sits.
//!
//! **A search on a bound needs a polish.** The bound is loose on a plateau of
//! near-optimal clearance, so the incumbent is the best *sampled* centre and can
//! sit a percent short of the true peak. A climb from the incumbent, then a
//! coarse rescan of the region with a climb from each strong sample, recovers it.
//! The rescan is what stops one walk from settling in a single basin.
use std::collections::BinaryHeap;

use super::super::boundary::Region;
use super::super::index::NearestObstacleIndex;
use super::super::input::ObstacleSet;
use crate::tuning;

const SQRT2: f64 = std::f64::consts::SQRT_2;

#[derive(Debug, Clone, Copy)]
struct GridCell {
    x: f64,
    y: f64,
    h_size: f64,
    max_dist: f64,
    /// Primitive ordinal achieving the nearest distance at this cell's centre,
    /// used to seed the next level so it starts from a near hit.
    support: usize,
}

impl PartialEq for GridCell {
    fn eq(&self, other: &Self) -> bool {
        self.max_dist == other.max_dist
    }
}

impl Eq for GridCell {}

impl PartialOrd for GridCell {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for GridCell {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        self.max_dist
            .partial_cmp(&other.max_dist)
            .unwrap_or(std::cmp::Ordering::Equal)
    }
}

/// One grid answer, before it is dressed as a public result.
#[derive(Debug, Clone, Copy)]
pub struct GridSolution {
    /// Centre x.
    pub x: f64,
    /// Centre y.
    pub y: f64,
    /// Squared clearance at the centre.
    pub radius_sq: f64,
    /// Cells the search popped.
    pub iterations: usize,
    /// True when the budget, not the bound, ended the search.
    pub exhausted_budget: bool,
}

/// Solve the largest empty circle by branch and bound over the region.
///
/// The returned radius is measured against the obstacle set at the returned
/// centre rather than read off a cell bound.
#[cfg_attr(feature = "hotpath", hotpath::measure)]
pub fn solve_grid(
    set: &ObstacleSet,
    nb: &NearestObstacleIndex,
    region: &Region,
    tolerance: f64,
) -> Option<GridSolution> {
    let (min_x, min_y, max_x, max_y) = region.bounds();

    let width = max_x - min_x;
    let height = max_y - min_y;
    let cell_size = width.max(height);

    if cell_size <= 0.0 || !tolerance.is_finite() || tolerance <= 0.0 {
        let (cx, cy) = centroid_of(set)?;
        if region.contains_closed(cx, cy) {
            let r_sq = set.clearance_sq(cx, cy);
            if r_sq.is_finite() && r_sq > 0.0 {
                return Some(GridSolution {
                    x: cx,
                    y: cy,
                    radius_sq: r_sq,
                    iterations: 0,
                    exhausted_budget: false,
                });
            }
        }
        return None;
    }

    // Half-size of a first-level cell: the box is covered by at most GRID_SIDE
    // cells along its longer axis and by however many are needed along the
    // shorter one, so a thin region is still covered end to end.
    let side = tuning::MEC_GRID_SIDE;
    let h_size = cell_size / (2.0 * side as f64);
    let step = 2.0 * h_size;
    let cells_x = ((width / step).ceil() as usize).max(1);
    let cells_y = ((height / step).ceil() as usize).max(1);

    let mut queue: BinaryHeap<GridCell> = BinaryHeap::with_capacity(cells_x * cells_y * 2);
    let mut best_x = f64::NAN;
    let mut best_y = f64::NAN;
    let mut best_dist = f64::NEG_INFINITY;

    for i in 0..cells_x {
        let cx = min_x + h_size + i as f64 * step;
        for j in 0..cells_y {
            let cy = min_y + h_size + j as f64 * step;
            let d = nearest_distance(nb, cx, cy);
            if region.contains_closed(cx, cy) && d > best_dist {
                best_x = cx;
                best_y = cy;
                best_dist = d;
            }
            queue.push(GridCell {
                x: cx,
                y: cy,
                h_size,
                max_dist: d + h_size * SQRT2,
                support: usize::MAX,
            });
        }
    }

    let mut iterations = 0usize;
    let mut exhausted_budget = false;

    while let Some(cell) = queue.pop() {
        iterations += 1;

        // Ordered by the bound, so once the best remaining cell cannot beat the
        // incumbent by more than the prune gap, no cell can. Runs before the
        // budget so a converged search ends on its own bound.
        if best_dist.is_finite() && cell.max_dist <= best_dist + tolerance {
            break;
        }

        if iterations > tuning::MEC_GRID_CELL_BUDGET {
            exhausted_budget = true;
            break;
        }

        let half = cell.h_size * 0.5;

        // Refining below the tolerance gains nothing: the child bound differs
        // from the parent's by at most `half * sqrt(2)`, already under the gap
        // the search resolves.
        if half * SQRT2 <= tolerance {
            continue;
        }

        let children = [
            (cell.x - half, cell.y - half),
            (cell.x + half, cell.y - half),
            (cell.x - half, cell.y + half),
            (cell.x + half, cell.y + half),
        ];

        // A child is pruned when `d + half*sqrt(2) <= best + tolerance`, which
        // rearranges to `d <= best - max(0, half*sqrt(2) - tolerance)`. Asking
        // the index for "the nearest distance, or a value under this bar" lets
        // the descent stop at the first primitive inside the bar instead of
        // proving the true minimum.
        let slack = half * SQRT2 - tolerance;
        let threshold = if best_dist.is_finite() && slack > 0.0 {
            Some(best_dist - slack)
        } else {
            None
        };

        for (cx, cy) in children {
            if cx + half < min_x || cx - half > max_x || cy + half < min_y || cy - half > max_y {
                continue;
            }

            let inside = region.contains_closed(cx, cy);
            let (d, support) = match threshold {
                Some(t) if inside => {
                    match nb.nearest_distance_sq_above(
                        cx,
                        cy,
                        (t.max(0.0)) * (t.max(0.0)),
                        Some(cell.support),
                    ) {
                        Some(v) => v,
                        // Proven within the bar: this child cannot improve.
                        None => continue,
                    }
                }
                _ => {
                    let Some(v) = nb.nearest_distance_sq(cx, cy) else {
                        continue;
                    };
                    v
                }
            };
            let raw = d.max(0.0).sqrt();
            let max_dist = raw + half * SQRT2;
            if best_dist.is_finite() && max_dist <= best_dist + tolerance {
                continue;
            }
            if inside && raw > best_dist {
                best_x = cx;
                best_y = cy;
                best_dist = raw;
            }
            queue.push(GridCell {
                x: cx,
                y: cy,
                h_size: half,
                max_dist,
                support,
            });
        }
    }

    if !best_dist.is_finite() || !best_x.is_finite() {
        return None;
    }

    // Polish. The search ranks cells by a bound, so its incumbent is the best
    // sampled centre; on a plateau that can sit short of the peak. The climb
    // fixes the local miss, the rescan fixes the basin choice.
    // The cheap climb always runs: it costs a few dozen probes and fixes the
    // local miss at no measurable price.
    let mut best = (best_x, best_y, best_dist);
    best = polish(nb, region, best, tolerance * 1e-4, h_size * 0.5);

    // The global rescan is the expensive recovery, so it only runs when the
    // search stopped on the budget. A search that converged on its own bound is
    // already within tolerance of the peak, and rescanning would add a thousand
    // distance queries to a case that took a handful.
    if exhausted_budget {
        best = polish_global(nb, region, best);
    }

    if !region.contains_closed(best.0, best.1) {
        return None;
    }

    let radius_sq = set.clearance_sq(best.0, best.1);
    if !radius_sq.is_finite() {
        return None;
    }

    Some(GridSolution {
        x: best.0,
        y: best.1,
        radius_sq: radius_sq.max(0.0),
        iterations,
        exhausted_budget,
    })
}

#[inline]
fn nearest_distance(nb: &NearestObstacleIndex, x: f64, y: f64) -> f64 {
    match nb.nearest_distance_sq(x, y) {
        Some((d, _)) if d.is_finite() && d > 0.0 => d.sqrt(),
        _ => 0.0,
    }
}

/// Climb the clearance field from a start point.
///
/// Sixteen directions rather than eight, because a ridge whose direction is not
/// axis- or diagonal-aligned is invisible to a compass walk. A successful round
/// grows the step, so a flat basin is crossed with long strides instead of
/// crawled; a failed one halves it. Only admissible, strictly better points are
/// accepted, so the walk cannot leave the region.
fn polish(
    nb: &NearestObstacleIndex,
    region: &Region,
    start: (f64, f64, f64),
    tol: f64,
    initial_step: f64,
) -> (f64, f64, f64) {
    // Sixteen directions at 22.5 degree steps, built at run time so the table
    // carries no rounded constants.
    let dirs: [(f64, f64); 16] = std::array::from_fn(|k| {
        let a = std::f64::consts::TAU * k as f64 / 16.0;
        (a.cos(), a.sin())
    });

    let (mut bx, mut by, mut bd) = start;
    if !region.contains_closed(bx, by) {
        return start;
    }

    let max_step = initial_step * 64.0;
    let mut step = initial_step;

    for _ in 0..2048 {
        if step <= tol || !step.is_finite() {
            break;
        }
        let mut moved = false;
        for pass in 0..2 {
            // Second pass rotates the probes by half the angular resolution.
            // Sixteen fixed directions leave an 11.25 degree blind gap, and a
            // crest lying inside that gap looks locally maximal to the first
            // pass; the rotated pass sees it. This is the difference between a
            // climb that converges to the peak and one that stops next to it.
            let phase = if pass == 0 {
                0.0
            } else {
                std::f64::consts::TAU / 32.0
            };
            for (dx, dy) in dirs {
                let (s, c) = phase.sin_cos();
                let rx = dx * c - dy * s;
                let ry = dx * s + dy * c;
                let nx = bx + rx * step;
                let ny = by + ry * step;
                if !region.contains_closed(nx, ny) {
                    continue;
                }
                let d = nearest_distance(nb, nx, ny);
                if d > bd {
                    bx = nx;
                    by = ny;
                    bd = d;
                    moved = true;
                }
            }
            if moved {
                break;
            }
        }
        if moved {
            step = (step * 1.5).min(max_step);
        } else {
            step *= 0.5;
        }
    }
    (bx, by, bd)
}

/// Coarse rescan of the region, then a climb from each strong sample.
///
/// The search ranks cells by a bound, so it can leave a whole basin unsampled,
/// and one walk from the incumbent then settles in whichever basin it began in.
/// Scanning coarsely, keeping the strongest samples, and climbing from each
/// covers every basin for a bounded cost. The scan supplies basins; the climb
/// supplies accuracy.
fn polish_global(
    nb: &NearestObstacleIndex,
    region: &Region,
    start: (f64, f64, f64),
) -> (f64, f64, f64) {
    let (min_x, min_y, max_x, max_y) = region.bounds();
    let span = (max_x - min_x).max(max_y - min_y);
    if !span.is_finite() || span <= 0.0 {
        return start;
    }

    let side = tuning::MEC_POLISH_SCAN_SIDE;
    let step_x = (max_x - min_x) / side as f64;
    let step_y = (max_y - min_y) / side as f64;
    if step_x <= 0.0 || step_y <= 0.0 {
        return start;
    }

    let mut seeds: Vec<(f64, f64, f64)> = Vec::with_capacity(tuning::MEC_POLISH_SEEDS);
    for i in 0..=side {
        for j in 0..=side {
            let x = min_x + step_x * i as f64;
            let y = min_y + step_y * j as f64;
            if !region.contains_closed(x, y) {
                continue;
            }
            let d = nearest_distance(nb, x, y);
            if d <= start.2 {
                continue;
            }
            seeds.push((x, y, d));
        }
    }
    seeds.sort_by(|a, b| b.2.partial_cmp(&a.2).unwrap_or(std::cmp::Ordering::Equal));
    seeds.truncate(tuning::MEC_POLISH_SEEDS);

    let mut best = start;
    for (x, y, d) in seeds {
        let (px, py, pd) = polish(nb, region, (x, y, d), span * 1e-12, step_x.max(step_y));
        if pd > best.2 {
            best = (px, py, pd);
        }
    }
    best
}

/// Centroid of the obstacle set: the mean of every point and segment midpoint.
fn centroid_of(set: &ObstacleSet) -> Option<(f64, f64)> {
    let mut sx = 0.0;
    let mut sy = 0.0;
    let mut n = 0usize;
    for p in &set.points {
        sx += p[0];
        sy += p[1];
        n += 1;
    }
    for i in 0..set.segments.len() {
        sx += (set.segments.ax[i] + set.segments.bx[i]) * 0.5;
        sy += (set.segments.ay[i] + set.segments.by[i]) * 0.5;
        n += 1;
    }
    (n > 0).then(|| (sx / n as f64, sy / n as f64))
}
