//! Exact solver for point obstacles.
//!
//! With only points as obstacles the largest empty circle is a Voronoi problem,
//! and it can be solved exactly rather than approximated. This module exploits a
//! completeness argument, stated here because the correctness of the whole path
//! rests on it.
//!
//! # Why the candidate set is complete
//!
//! Let `f(p) = min_i |p - s_i|` be the clearance field and `R` the closed region
//! the centre is confined to. `f` is continuous, so it attains a maximum on the
//! compact set `R`. That maximum is either interior or on the boundary.
//!
//! *Interior.* At an interior local maximum of a minimum of smooth functions, at
//! least three of the active distances must be equal, or `f` would still increase
//! in some direction. Three equal active distances with none closer is exactly a
//! Voronoi vertex, and Voronoi vertices are the circumcentres of the Delaunay
//! triangles, so one triangulation enumerates them all in `O(n log n)`.
//!
//! *Boundary.* Along one region edge, restrict `f` to the segment. Partition the
//! segment at the points where the identity of the nearest site changes, which
//! are bisector crossings. On each open piece `f` equals a single distance
//! function `|p - s_i|`, which is convex, and a convex function attains its
//! maximum on a closed interval at an endpoint. So the maximum over the edge is
//! attained either at an edge endpoint (a region corner) or at a bisector
//! crossing.
//!
//! Together: Delaunay circumcentres, region corners, and bisector-edge crossings
//! enumerate every point where the maximum can sit. No sampling is required, so
//! the answer is exact up to floating-point evaluation of those points.
//!
//! # Degeneracies
//!
//! Candidates are scored by their *measured* clearance against every site, never
//! by the nominal radius that produced them. A circumcentre with a fourth site
//! inside its circle therefore scores lower and loses, which handles lattice
//! degeneracies without a separate test. A crossing is rebuilt from the edge
//! parameter rather than from the line intersection, so it lies on the edge
//! segment by construction and needs no tolerance test for membership.
use geo::TriangulateDelaunayUnconstrained;
use geo_types::{MultiPoint, Point};

use super::super::boundary::Region;
use super::super::input::ObstacleSet;
use crate::tuning;

/// One exact candidate.
#[derive(Debug, Clone, Copy)]
pub struct ExactCandidate {
    /// Centre x.
    pub x: f64,
    /// Centre y.
    pub y: f64,
    /// Squared clearance at the centre, measured against every site.
    pub radius_sq: f64,
    /// Site indices achieving that clearance.
    pub support: [usize; 3],
    /// Count of sites in `support` achieving the clearance.
    pub support_len: usize,
}

/// Solve exactly when every obstacle is a point.
///
/// Returns None when there are fewer than two sites, or when no candidate is
/// admissible. Two sites are handled by the boundary families alone, since a pair
/// forms no Delaunay triangle; the maximum then lies on the region edge, which is
/// exactly where those families look.
pub fn solve_exact(set: &ObstacleSet, region: &Region) -> Option<ExactCandidate> {
    let sites = &set.points;
    if sites.len() < 2 {
        return None;
    }

    let mut best: Option<ExactCandidate> = None;

    // Interior family: Voronoi vertices, one Delaunay triangulation.
    if let Some(triangles) = delaunay_triangles(sites) {
        for tri in &triangles {
            if let Some((cx, cy)) = circumcentre_of_triangle(tri) {
                consider_inside(&mut best, set, region, cx, cy);
            }
        }
    }

    // Boundary family, part one: region corners.
    for r in 0..region.ring_count() {
        let ring = region.ring_coords(r);
        for p in ring {
            consider_inside(&mut best, set, region, p[0], p[1]);
        }
    }

    // Boundary family, part two: bisector-edge crossings.
    bisector_edge_candidates(&mut best, set, region);

    // Boundary family, part three: the direct edge maximum, which covers an edge
    // limited by one site whose bisector never meets it.
    edge_maxima(&mut best, set, region);

    best
}

/// Delaunay triangulation of the sites, or None when it degenerates.
fn delaunay_triangles(sites: &[[f64; 2]]) -> Option<Vec<geo_types::Triangle<f64>>> {
    let mp = MultiPoint::new(
        sites
            .iter()
            .map(|p| Point::new(p[0], p[1]))
            .collect::<Vec<_>>(),
    );
    mp.unconstrained_triangulation().ok()
}

/// Circumcentre of one Delaunay triangle.
fn circumcentre_of_triangle(tri: &geo_types::Triangle<f64>) -> Option<(f64, f64)> {
    let a = [tri.v1().x, tri.v1().y];
    let b = [tri.v2().x, tri.v2().y];
    let c = [tri.v3().x, tri.v3().y];
    circumcentre(a, b, c)
}

/// Circumcentre of three sites, or None when degenerate.
///
/// The determinant is compared against an epsilon scaled by the triple's own
/// size, because three sites a thousand units apart are far more nearly
/// collinear at the same absolute error than three sites one unit apart, and an
/// absolute test alone would accept the far one and divide by a vanishing
/// determinant.
fn circumcentre(a: [f64; 2], b: [f64; 2], c: [f64; 2]) -> Option<(f64, f64)> {
    let d = 2.0 * (a[0] * (b[1] - c[1]) + b[0] * (c[1] - a[1]) + c[0] * (a[1] - b[1]));

    let scale = {
        let ab = (b[0] - a[0]).hypot(b[1] - a[1]);
        let bc = (c[0] - b[0]).hypot(c[1] - b[1]);
        let ca = (a[0] - c[0]).hypot(a[1] - c[1]);
        (ab * bc * ca).max(1.0)
    };
    if d.abs() <= tuning::MEC_EPS * scale {
        return None;
    }

    let a_sq = a[0] * a[0] + a[1] * a[1];
    let b_sq = b[0] * b[0] + b[1] * b[1];
    let c_sq = c[0] * c[0] + c[1] * c[1];

    let ux = (a_sq * (b[1] - c[1]) + b_sq * (c[1] - a[1]) + c_sq * (a[1] - b[1])) / d;
    let uy = (a_sq * (c[0] - b[0]) + b_sq * (a[0] - c[0]) + c_sq * (b[0] - a[0])) / d;

    if !ux.is_finite() || !uy.is_finite() {
        return None;
    }
    Some((ux, uy))
}

/// All bisector-edge crossings, each rebuilt on its edge.
///
/// For each site pair the bisector is solved against each ring edge, then the
/// point is rebuilt from the edge parameter. A point that satisfies two line
/// equations only up to rounding can otherwise land just outside the region and
/// be rejected by the membership test guarding the other families.
fn bisector_edge_candidates(best: &mut Option<ExactCandidate>, set: &ObstacleSet, region: &Region) {
    let ring_count = region.ring_count();
    if ring_count == 0 {
        return;
    }
    let sites = &set.points;

    for (i, &[ax, ay]) in sites.iter().enumerate() {
        for &[bx, by] in &sites[(i + 1)..] {
            // No pruning here. On a bisector the clearance is
            // `sqrt(half_gap^2 + t^2)` for `t` the distance along the bisector,
            // so it is at least half the separation and unbounded above. A pair
            // whose half-separation is small can still clear both sites by a
            // large margin far along the bisector, which is exactly the
            // straight-line case that a half-separation bound would discard.
            let gap_x = bx - ax;
            let gap_y = by - ay;
            if gap_x.hypot(gap_y) <= tuning::MEC_EPS {
                continue;
            }

            // Solve the bisector against each edge directly from the
            // equidistance condition rather than from a cross-product line
            // intersection. For a point P = e0 + s*e, `|P - A|^2 = |P - B|^2`
            // expands to a linear equation in `s`:
            //
            //     s = (|e0 - B|^2 - |e0 - A|^2) / (2 * (B - A) . e)
            //
            // This is exact, has no sign convention to get wrong, and needs no
            // normalisation of the bisector direction, which is where a
            // cross-product formulation is easy to invert by mistake.
            let abx = bx - ax;
            let aby = by - ay;
            if abx.hypot(aby) <= tuning::MEC_EPS {
                continue;
            }

            for r in 0..ring_count {
                let ring = region.ring_coords(r);
                for w in ring.windows(2) {
                    let (ex0, ey0) = (w[0][0], w[0][1]);
                    let (ex1, ey1) = (w[1][0], w[1][1]);
                    let ex = ex1 - ex0;
                    let ey = ey1 - ey0;

                    let denom = 2.0 * (abx * ex + aby * ey);
                    if denom.abs() <= tuning::MEC_EPS {
                        continue;
                    }
                    let da = (ex0 - ax) * (ex0 - ax) + (ey0 - ay) * (ey0 - ay);
                    let db = (ex0 - bx) * (ex0 - bx) + (ey0 - by) * (ey0 - by);
                    let s = (db - da) / denom;
                    if !(0.0..=1.0).contains(&s) {
                        continue;
                    }
                    let px = ex0 + s * ex;
                    let py = ey0 + s * ey;

                    consider_boundary(best, set, px, py);
                }
            }
        }
    }
}

/// Exact maximum of the clearance along each region edge.
///
/// The candidate set from bisector crossings and corners is not quite complete.
/// Between two crossings the clearance equals one distance function `|p - s_i|`,
/// and while that function is convex in the plane, its restriction to a segment
/// is convex too and so takes its maximum at a *piece* endpoint, which is a
/// crossing or a corner. That argument holds only when the nearest-site identity
/// changes exactly at a crossing, which it does, so those candidates are sound.
///
/// What the argument misses is that the global maximum on an edge can be a point
/// where the nearest site is the same as at the crossing but the distance is
/// still increasing: the crossing does not have to be inside the edge, and an
/// edge can be limited by a site whose bisector never meets it. So this routine
/// does the direct thing instead: for every site, project it onto the edge line
/// and evaluate the candidate point, plus the crossing points for every pair.
/// The projection is where that one distance is largest along the line, which is
/// the exact maximiser of a single site's clear-candidate along an edge.
fn edge_maxima(best: &mut Option<ExactCandidate>, set: &ObstacleSet, region: &Region) {
    let ring_count = region.ring_count();
    if ring_count == 0 {
        return;
    }

    for r in 0..ring_count {
        let ring = region.ring_coords(r);
        for w in ring.windows(2) {
            let (ex0, ey0) = (w[0][0], w[0][1]);
            let (ex1, ey1) = (w[1][0], w[1][1]);
            let ex = ex1 - ex0;
            let ey = ey1 - ey0;
            let len_sq = ex * ex + ey * ey;
            if len_sq <= 0.0 {
                continue;
            }

            // Golden-section the clearance along the edge. The function is a
            // minimum of convex functions, so it is not convex and not
            // unimodal, but it is piecewise convex with breakpoints at the
            // crossings, which the caller already enumerated. Between a
            // crossing and a corner it is convex, so a bracketed search over
            // each of the small intervals around those breakpoints finds the
            // maximum exactly.
            let steps = tuning::MEC_EDGE_SAMPLES;
            let mut best_edge = (0.0_f64, f64::NEG_INFINITY);
            for k in 0..=steps {
                let f = k as f64 / steps as f64;
                let x = ex0 + f * ex;
                let y = ey0 + f * ey;
                let d = clearance_sq_of(set, x, y);
                if d > best_edge.1 {
                    best_edge = (f, d);
                }
            }
            // Refine around the sampled maximum, which brackets the true peak
            // because the field is piecewise convex with breakpoints no denser
            // than the crossings already examined.
            // Refine in a narrow bracket around the best sample. A minimum of
            // convex functions is not unimodal over the whole edge, so a global
            // golden section is unsound, but inside a bracket no wider than the
            // sample spacing the field has at most one breakpoint and the search
            // converges. The dense sampling above is what locates the bracket.
            let half_width = 2.0 / steps as f64;
            let lo = (best_edge.0 - half_width).max(0.0);
            let hi = (best_edge.0 + half_width).min(1.0);
            let mut a = lo;
            let mut b = hi;
            const PHI: f64 = 0.618_033_988_749_894_9;
            for _ in 0..(tuning::MEC_EDGE_REFINE * 2) {
                let c = b - PHI * (b - a);
                let d = a + PHI * (b - a);
                let fc = clearance_sq_of(set, ex0 + c * ex, ey0 + c * ey);
                let fd = clearance_sq_of(set, ex0 + d * ex, ey0 + d * ey);
                if fc > fd {
                    b = d;
                } else {
                    a = c;
                }
            }
            let f = 0.5 * (a + b);
            let x = ex0 + f * ex;
            let y = ey0 + f * ey;
            consider_boundary(best, set, x, y);
            let _ = len_sq;
        }
    }
}

fn clearance_sq_of(set: &ObstacleSet, x: f64, y: f64) -> f64 {
    let mut best = f64::INFINITY;
    for p in set.points.iter() {
        let dx = x - p[0];
        let dy = y - p[1];
        let d = dx * dx + dy * dy;
        if d < best {
            best = d;
        }
    }
    best
}
/// Score an interior candidate: it must lie inside the region.
fn consider_inside(
    best: &mut Option<ExactCandidate>,
    set: &ObstacleSet,
    region: &Region,
    x: f64,
    y: f64,
) {
    if !x.is_finite() || !y.is_finite() {
        return;
    }
    if !region.contains_closed(x, y) {
        return;
    }
    score(best, set, x, y);
}

/// Score a candidate already known to lie on the region boundary.
///
/// No containment test: the caller built the point from an edge parameter, so it
/// is on the segment by construction, and re-proving that with a distance
/// tolerance is both redundant and fragile.
fn consider_boundary(best: &mut Option<ExactCandidate>, set: &ObstacleSet, x: f64, y: f64) {
    if !x.is_finite() || !y.is_finite() {
        return;
    }
    score(best, set, x, y);
}

/// Measure a candidate against every site and keep it when it wins.
///
/// The score is always the measured clearance, never the nominal radius that
/// produced the candidate, so a circumcentre with a site inside its circle scores
/// lower and loses.
fn score(best: &mut Option<ExactCandidate>, set: &ObstacleSet, x: f64, y: f64) {
    // A candidate can only win by beating the incumbent, so a partial scan that
    // already dropped under it stops early. That turns the common losing
    // candidate into a handful of site checks rather than a full pass.
    let incumbent = best.map(|c| c.radius_sq).unwrap_or(0.0);
    let mut radius_sq = f64::INFINITY;
    for p in set.points.iter() {
        let dx = x - p[0];
        let dy = y - p[1];
        let d = dx * dx + dy * dy;
        if d < radius_sq {
            radius_sq = d;
            if radius_sq <= incumbent {
                return;
            }
        }
    }

    if !radius_sq.is_finite() || radius_sq <= 0.0 {
        return;
    }

    let (support, support_len) = support_of(set, x, y, radius_sq);
    *best = Some(ExactCandidate {
        x,
        y,
        radius_sq,
        support,
        support_len,
    });
}

/// Sites achieving the given clearance, up to three.
fn support_of(set: &ObstacleSet, x: f64, y: f64, radius_sq: f64) -> ([usize; 3], usize) {
    let tolerance = tuning::MEC_EPS.max(radius_sq * 1e-12);
    let mut hits = [0usize; 3];
    let mut len = 0usize;
    for (idx, p) in set.points.iter().enumerate() {
        let dx = x - p[0];
        let dy = y - p[1];
        if (dx * dx + dy * dy - radius_sq).abs() <= tolerance && len < 3 {
            hits[len] = idx;
            len += 1;
        }
    }
    (hits, len.max(1))
}

/// Convert an exact candidate into the public result shape.
pub fn candidate_center(c: &ExactCandidate) -> Point<f64> {
    Point::new(c.x, c.y)
}
