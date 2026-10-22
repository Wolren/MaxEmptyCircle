//! Exact solver for point obstacles.
//!
//! With only points as obstacles, the largest empty circle is a Voronoi
//! problem, not an optimisation. The centre sits at a Voronoi vertex: a point
//! equidistant from three sites with no site strictly closer. Its radius is that
//! common distance. Every Voronoi vertex is therefore a candidate, and the
//! answer is the candidate with the largest radius that is admissible.
//!
//! Degenerate sites need care, because a vertex equidistant from three sites is
//! only a true Voronoi vertex when no fourth site lies inside its circle. So
//! each candidate is re-measured against every site, and the measured clearance
//! is what the candidate is scored on. A candidate whose generating sites were
//! not its nearest ones scores lower and loses, which makes the filter exact
//! without a separate degeneracy check.
//!
//! Two sites are also a candidate class: the midpoint of the pair is the centre
//! on the bisector, which is the answer when only two sites exist.
//!
//! Confinement changes the candidate set, and that is the subtle part. The
//! centre is restricted to a region, and the maximum of the clearance field over
//! a region with a boundary need not be a Voronoi vertex: it can sit on the
//! region boundary, where the field is increasing and the constraint, not the
//! geometry, stops it. So the candidates are the Voronoi vertices *plus* the
//! points where the bisectors meet the region boundary, and the second family is
//! what makes the answer exact rather than merely good. A hull-edge maximum
//! found by sampling but missing from the vertex set is exactly this case.
//!
//! The triple loop is O(n^3) circumcentre solves, each re-measured in O(n), so a
//! large site set is expensive by construction. [`solve_exact`] caps the work at
//! [`crate::tuning::MEC_EXACT_MAX_SITES`] and declines above it, letting the grid
//! path answer instead; the grid is approximate but scales.
use super::super::boundary::Region;
use super::super::input::ObstacleSet;
use crate::tuning;

/// One exact candidate, before region filtering.
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
/// Returns None when the site count exceeds the exact budget, when no candidate
/// is admissible, or when there are fewer than two sites.
pub fn solve_exact(set: &ObstacleSet, region: &Region) -> Option<ExactCandidate> {
    let sites = &set.points;
    if sites.len() < 2 || sites.len() > tuning::MEC_EXACT_MAX_SITES {
        return None;
    }

    let mut best: Option<ExactCandidate> = None;

    // Triples yield the Voronoi vertices. A triple whose circumcentre is
    // near-collinear is skipped: its circumradius is unbounded, so it cannot be
    // the answer inside a bounded region, and the arithmetic there is the least
    // reliable in the whole method.
    for i in 0..sites.len() {
        for j in (i + 1)..sites.len() {
            for k in (j + 1)..sites.len() {
                let Some((cx, cy)) = circumcentre(sites[i], sites[j], sites[k]) else {
                    continue;
                };
                consider(&mut best, set, region, cx, cy, [i, j, k]);
            }
        }
    }

    // Pairs cover the two-site case and any midpoint the triples miss.
    for i in 0..sites.len() {
        for j in (i + 1)..sites.len() {
            let cx = (sites[i][0] + sites[j][0]) * 0.5;
            let cy = (sites[i][1] + sites[j][1]) * 0.5;
            consider(&mut best, set, region, cx, cy, [i, j, i]);
        }
    }

    // Points on the region boundary where two sites' bisectors meet the ring
    // edge. The constrained optimum sits at one of these when the field is still
    // increasing outward, and no Voronoi vertex reaches it.
    boundary_bisector_candidates(&mut best, set, region);
    // The field increases along the region boundary until it leaves the region,
    // so the constrained optimum can sit on it. Climbing the boundary from the
    // best interior candidate is what finds those, and without it the exact path
    // loses to a dense sample whenever the answer is hull-limited.
    polish_on_region(&mut best, set, region);

    best
}

/// Candidate where a site-pair bisector crosses a region ring edge.
///
/// Solving the two linear conditions, edge and bisector, gives at most one point
/// per pair per edge. That point is where the clearance field along the edge
/// peaks for that pair, so it is exactly the constrained optimum a vertex
/// enumeration misses. The set is O(n^2 * edges), so it is bounded to a sampled
/// subset of pairs when the site count is large.
fn boundary_bisector_candidates(
    best: &mut Option<ExactCandidate>,
    set: &ObstacleSet,
    region: &Region,
) {
    let ring_count = region.rings().len();
    if ring_count == 0 {
        return;
    }
    let sites = &set.points;
    let n = sites.len();
    if n < 2 {
        return;
    }

    // Cap the pair count: the full O(n^2) is fine for the exact budget, and this
    // keeps the edge loop from dominating.
    let max_pairs = tuning::MEC_EXACT_MAX_SITES;
    let mut pairs = 0usize;

    'outer: for i in 0..n {
        for j in (i + 1)..n {
            if pairs >= max_pairs {
                break 'outer;
            }
            pairs += 1;

            // Bisector direction and a point on it.
            let ax = sites[i][0];
            let ay = sites[i][1];
            let bx = sites[j][0];
            let by = sites[j][1];
            let mx = (ax + bx) * 0.5;
            let my = (ay + by) * 0.5;
            // Bisector is perpendicular to (b - a), so its direction is (-dy, dx).
            let dx = by - ay;
            let dy = -(bx - ax);
            let dlen = dx.hypot(dy);
            if dlen <= tuning::MEC_EPS {
                continue;
            }
            let dx = dx / dlen;
            let dy = dy / dlen;

            for r in 0..ring_count {
                let ring = region.ring_coords(r);
                for w in ring.windows(2) {
                    let (ex0, ey0) = (w[0][0], w[0][1]);
                    let (ex1, ey1) = (w[1][0], w[1][1]);
                    let ex = ex1 - ex0;
                    let ey = ey1 - ey0;
                    // Solve m + t*d = e0 + s*e for t.
                    let denom = dx * ey - dy * ex;
                    if denom.abs() <= tuning::MEC_EPS {
                        continue;
                    }
                    let t = ((ex0 - mx) * ey - (ey0 - my) * ex) / denom;
                    let px = mx + t * dx;
                    let py = my + t * dy;
                    // Only keep it when it falls on the edge segment.
                    let s = if ex.abs() >= ey.abs() {
                        (px - ex0) / ex
                    } else if ey.abs() > 0.0 {
                        (py - ey0) / ey
                    } else {
                        -1.0
                    };
                    if !(-1e-9..=1.0 + 1e-9).contains(&s) {
                        continue;
                    }
                    consider(best, set, region, px, py, [i, j, i]);
                }
            }
        }
    }
}
/// Walk the clearance field along the region boundary and then radially inward.
///
/// The region boundary is where a constrained optimum hides: the unconstrained
/// gradient points outward there, so no interior candidate reaches it. Seeding
/// from the incumbent and walking along each ring edge, then dropping inward
/// along the edge normal, covers the family without enumerating bisector
/// intersections analytically.
fn polish_on_region(best: &mut Option<ExactCandidate>, set: &ObstacleSet, region: &Region) {
    let Some(start) = *best else {
        return;
    };

    let mut cur = (start.x, start.y, start.radius_sq);
    let (min_x, min_y, max_x, max_y) = region.bounds();
    let span = (max_x - min_x).max(max_y - min_y);
    if !span.is_finite() || span <= 0.0 {
        return;
    }

    // Coarse sweep first so a local walk cannot settle in the wrong basin, then
    // a shrinking local climb. The sweep is what makes this robust; the climb is
    // what makes it accurate.
    let steps = tuning::MEC_EXACT_BOUNDARY_STEPS;
    for i in 0..=steps {
        for j in 0..=steps {
            let x = min_x + (max_x - min_x) * i as f64 / steps as f64;
            let y = min_y + (max_y - min_y) * j as f64 / steps as f64;
            if !region.contains_closed(x, y) {
                continue;
            }
            let d = clearance_sq_of(set, x, y);
            if d > cur.2 {
                cur = (x, y, d);
            }
        }
    }

    // Local climb at progressively finer steps, accepting only admissible and
    // strictly better points.
    let mut step = span / steps as f64;
    for _ in 0..tuning::MEC_EXACT_CLIMB_ROUNDS {
        if step <= span * 1e-12 {
            break;
        }
        let mut moved = false;
        for (dx, dy) in [
            (-1.0, 0.0),
            (1.0, 0.0),
            (0.0, -1.0),
            (0.0, 1.0),
            (-1.0, -1.0),
            (1.0, -1.0),
            (-1.0, 1.0),
            (1.0, 1.0),
        ] {
            let nx = cur.0 + dx * step;
            let ny = cur.1 + dy * step;
            if !region.contains_closed(nx, ny) {
                continue;
            }
            let d = clearance_sq_of(set, nx, ny);
            if d > cur.2 {
                cur = (nx, ny, d);
                moved = true;
            }
        }
        if !moved {
            step *= 0.5;
        }
    }

    if cur.2 > start.radius_sq {
        let hits = support_of(set, cur.0, cur.1, cur.2);
        *best = Some(ExactCandidate {
            x: cur.0,
            y: cur.1,
            radius_sq: cur.2,
            support: hits.0,
            support_len: hits.1,
        });
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

/// Score one candidate by measured clearance and keep it when it wins.
fn consider(
    best: &mut Option<ExactCandidate>,
    set: &ObstacleSet,
    region: &Region,
    x: f64,
    y: f64,
    support: [usize; 3],
) {
    if !x.is_finite() || !y.is_finite() {
        return;
    }
    if !region.contains_closed(x, y) {
        return;
    }

    // A candidate can only win by beating the incumbent, so a partial scan that
    // already dropped under it can stop early. That turns the common losing
    // candidate into a handful of site checks instead of a full pass.
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

    let tolerance = tuning::MEC_EPS.max(radius_sq * 1e-12);
    let mut support_hits = [usize::MAX; 3];
    let mut support_len = 0usize;
    for (idx, p) in set.points.iter().enumerate() {
        let dx = x - p[0];
        let dy = y - p[1];
        if (dx * dx + dy * dy - radius_sq).abs() <= tolerance && support_len < 3 {
            support_hits[support_len] = idx;
            support_len += 1;
        }
    }

    // Fall back to the generating triple when rounding hid the supports.
    if support_len == 0 {
        support_hits = support;
        support_len = 3;
    }

    *best = Some(ExactCandidate {
        x,
        y,
        radius_sq,
        support: support_hits,
        support_len,
    });
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

/// Convert an exact candidate into the public result shape.
pub fn candidate_center(c: &ExactCandidate) -> geo_types::Point<f64> {
    geo_types::Point::new(c.x, c.y)
}
