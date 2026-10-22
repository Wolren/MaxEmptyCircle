//! Point-in-boundary index and the region a centre is confined to.
//!
//! Containment uses even-odd ray casting. Ring bounding boxes short-circuit the
//! rings that cannot contain the query, which matters on boundaries with holes:
//! the point only pays for the rings whose box it falls inside.
//!
//! [`Region`] is the well-posedness device. An empty circle between two finite
//! obstacles is unbounded: moving far enough away from both clears them by an
//! arbitrary distance. So the centre must be confined, and the region is what
//! confines it. A caller-supplied boundary is used when present; otherwise the
//! convex hull of the obstacles is used, which is the GEOS default.
use geo::ConvexHull;
use geo_types::{MultiPoint, Point};

use super::input::{BoundarySet, ObstacleSet};
use crate::MecError;

/// Point-in-boundary test with hole support.
#[derive(Debug, Clone)]
pub struct BoundaryIndex {
    rings: Vec<Ring>,
}

#[derive(Debug, Clone)]
pub(crate) struct Ring {
    coords: Vec<[f64; 2]>,
    bbox: (f64, f64, f64, f64),
}

impl Ring {
    fn new(coords: Vec<[f64; 2]>) -> Self {
        let mut bbox = (
            f64::INFINITY,
            f64::INFINITY,
            f64::NEG_INFINITY,
            f64::NEG_INFINITY,
        );
        for p in &coords {
            bbox.0 = bbox.0.min(p[0]);
            bbox.1 = bbox.1.min(p[1]);
            bbox.2 = bbox.2.max(p[0]);
            bbox.3 = bbox.3.max(p[1]);
        }
        Self { coords, bbox }
    }
}

impl BoundaryIndex {
    /// Build the containment index from a boundary set.
    pub fn new(set: &BoundarySet) -> Self {
        let rings = set
            .rings
            .iter()
            .map(|ring| Ring::new(ring.clone()))
            .collect();
        Self { rings }
    }

    /// Number of rings in the boundary.
    pub fn ring_count(&self) -> usize {
        self.rings.len()
    }

    /// Bounds of every ring as (min_x, min_y, max_x, max_y), or None.
    pub fn bounds(&self) -> Option<(f64, f64, f64, f64)> {
        let mut out: Option<(f64, f64, f64, f64)> = None;
        for ring in &self.rings {
            let (min_x, min_y, max_x, max_y) = ring.bbox;
            out = Some(match out {
                Some((a, b, c, d)) => (a.min(min_x), b.min(min_y), c.max(max_x), d.max(max_y)),
                None => (min_x, min_y, max_x, max_y),
            });
        }
        out
    }

    #[inline]
    /// True when (x, y) is inside the boundary, with hole interiors excluded.
    pub fn contains_xy(&self, x: f64, y: f64) -> bool {
        let mut inside = false;
        for ring in &self.rings {
            let (min_x, min_y, max_x, max_y) = ring.bbox;
            if x < min_x || x > max_x || y < min_y || y > max_y {
                continue;
            }
            if point_in_ring(x, y, &ring.coords) {
                inside = !inside;
            }
        }
        inside
    }
}

/// The closed region a circle centre is confined to.
///
/// A centre outside the region is not a candidate. Derived once per solve and
/// shared by both solve paths.
#[derive(Debug, Clone)]
pub struct Region {
    rings: Vec<Ring>,
    bounds: (f64, f64, f64, f64),
}

impl Region {
    /// Build the region: the caller boundary when present, else the obstacle
    /// convex hull.
    ///
    /// A degenerate hull, which is a single point or a set of collinear points,
    /// encloses no area and therefore confines nothing: the empty circle is
    /// unbounded and the solve is rejected with a message saying so.
    pub fn new(set: &ObstacleSet, boundary: Option<&BoundaryIndex>) -> Result<Self, MecError> {
        if let Some(b) = boundary {
            let bounds = b
                .bounds()
                .ok_or_else(|| MecError::InvalidInput("boundary has no coordinates".to_string()))?;
            return Ok(Self {
                rings: b.rings.clone(),
                bounds,
            });
        }

        let vertices = set.vertices();
        if vertices.len() < 3 {
            return Err(MecError::InvalidInput(
                "obstacles enclose no area: pass a polygonal boundary so the empty circle is bounded"
                    .to_string(),
            ));
        }

        let hull = MultiPoint::new(
            vertices
                .iter()
                .map(|p| Point::new(p[0], p[1]))
                .collect::<Vec<_>>(),
        )
        .convex_hull();

        let ring: Vec<[f64; 2]> = hull.exterior().0.iter().map(|c| [c.x, c.y]).collect();
        // A collinear hull comes back as a degenerate ring with no area, which
        // confines nothing, so the same rejection applies.
        if ring.len() < 4 || signed_area(&ring).abs() <= 1e-12 {
            return Err(MecError::InvalidInput(
                "obstacles are collinear: pass a polygonal boundary so the empty circle is bounded"
                    .to_string(),
            ));
        }

        let bounds = set
            .bounds()
            .ok_or_else(|| MecError::InvalidInput("obstacle set is empty".to_string()))?;

        Ok(Self {
            rings: vec![Ring::new(ring)],
            bounds,
        })
    }

    /// Search box the solver tiles: the region bounds.
    pub fn bounds(&self) -> (f64, f64, f64, f64) {
        self.bounds
    }

    /// The region rings, for a solver that walks the boundary.
    pub(crate) fn rings(&self) -> &[Ring] {
        &self.rings
    }

    /// True when the region has any ring, so walking the boundary is meaningful.
    pub fn has_rings(&self) -> bool {
        !self.rings.is_empty()
    }

    /// The ring coordinates, in order.
    pub fn ring_coords(&self, idx: usize) -> &[[f64; 2]] {
        &self.rings[idx].coords
    }

    /// True when (x, y) is inside the region or on its ring.
    ///
    /// The region is closed, so a point exactly on a ring is admissible. The
    /// plain containment test is an even-odd ray cast, which is ambiguous on the
    /// ring itself and rejects it, and that matters here: with obstacles, the
    /// constrained optimum of the clearance field sits exactly on the region
    /// boundary whenever the field is still increasing outward there.
    pub fn contains_closed(&self, x: f64, y: f64) -> bool {
        if self.contains(x, y) {
            return true;
        }
        let eps = region_eps(&self.bounds);
        for ring in &self.rings {
            for w in ring.coords.windows(2) {
                if point_segment_distance_sq(w[0], w[1], x, y) <= eps * eps {
                    return true;
                }
            }
        }
        false
    }

    #[inline]
    /// True when (x, y) lies inside the region.
    pub fn contains(&self, x: f64, y: f64) -> bool {
        let (min_x, min_y, max_x, max_y) = self.bounds;
        if x < min_x || x > max_x || y < min_y || y > max_y {
            return false;
        }
        let mut inside = false;
        for ring in &self.rings {
            if point_in_ring(x, y, &ring.coords) {
                inside = !inside;
            }
        }
        inside
    }
}

/// Signed area of a closed ring, by the shoelace formula.
fn signed_area(ring: &[[f64; 2]]) -> f64 {
    let n = ring.len();
    if n < 3 {
        return 0.0;
    }
    let mut sum = 0.0;
    for i in 0..n {
        let a = ring[i];
        let b = ring[(i + 1) % n];
        sum += a[0] * b[1] - b[0] * a[1];
    }
    sum * 0.5
}

fn point_in_ring(x: f64, y: f64, ring: &[[f64; 2]]) -> bool {
    #[cfg(not(feature = "shewchuk"))]
    {
        let mut inside = false;
        for w in ring.windows(2) {
            let (ax, ay) = (w[0][0], w[0][1]);
            let (bx, by) = (w[1][0], w[1][1]);
            if (ay > y) != (by > y) {
                let intersect_x = ax + (y - ay) * (bx - ax) / (by - ay);
                if x < intersect_x {
                    inside = !inside;
                }
            }
        }
        inside
    }
    #[cfg(feature = "shewchuk")]
    {
        use robust::{orient2d, Coord};
        let q = Coord { x, y };
        let mut inside = false;
        for w in ring.windows(2) {
            let a = Coord {
                x: w[0][0],
                y: w[0][1],
            };
            let b = Coord {
                x: w[1][0],
                y: w[1][1],
            };
            // Same half-open rule as the ray cast: a vertex on the ray counts once.
            if (a.y > y) != (b.y > y) {
                let o = orient2d(a, b, q);
                let left = if b.y > a.y { o > 0.0 } else { o < 0.0 };
                if left {
                    inside = !inside;
                }
            }
        }
        inside
    }
}

/// Tolerance for treating a point as lying on the region boundary.
fn region_eps(bounds: &(f64, f64, f64, f64)) -> f64 {
    let (min_x, min_y, max_x, max_y) = *bounds;
    let span = (max_x - min_x).max(max_y - min_y);
    if span.is_finite() && span > 0.0 {
        span * 1e-12
    } else {
        1e-12
    }
}

#[inline]
fn point_segment_distance_sq(a: [f64; 2], b: [f64; 2], x: f64, y: f64) -> f64 {
    let dx = b[0] - a[0];
    let dy = b[1] - a[1];
    let len_sq = dx * dx + dy * dy;
    let t = if len_sq > 0.0 {
        (((x - a[0]) * dx + (y - a[1]) * dy) / len_sq).clamp(0.0, 1.0)
    } else {
        0.0
    };
    let ex = x - (a[0] + t * dx);
    let ey = y - (a[1] + t * dy);
    ex * ex + ey * ey
}
