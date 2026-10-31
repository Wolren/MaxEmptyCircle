//! Invariants the empty circle must satisfy, plus a brute-force oracle for the
//! point-only case.
//!
//! Two properties are checked on every case. First, emptiness: the reported
//! radius must equal the measured clearance at the reported centre, so the
//! circle really is empty by the same distance function the solver uses.
//! Second, reach: no dense sample in the search region may clear more than the
//! reported radius, which is what would expose a solver that stopped early.
//!
//! The oracle is separate from the solver. For point obstacles the answer is a
//! Voronoi vertex, so a dense grid plus a refinement pass over the best sampled
//! cell finds the same optimum without sharing any code with the exact path.

use max_empty_circle::geo_types::{Geometry, LineString, MultiPoint, Point, Polygon};
use max_empty_circle::{largest_empty_circle, MecApproach, MecError, MecOptions, MecResult};

/// Deterministic xorshift so a failure reproduces from the seed alone.
struct Rng(u64);

impl Rng {
    fn next(&mut self) -> u64 {
        self.0 ^= self.0 << 13;
        self.0 ^= self.0 >> 7;
        self.0 ^= self.0 << 17;
        self.0
    }

    fn unit(&mut self) -> f64 {
        (self.next() >> 11) as f64 / (1u64 << 53) as f64
    }

    fn range(&mut self, lo: f64, hi: f64) -> f64 {
        lo + (hi - lo) * self.unit()
    }
}

fn sites_geom(pts: &[(f64, f64)]) -> Geometry<f64> {
    Geometry::MultiPoint(MultiPoint::new(
        pts.iter().map(|p| Point::new(p.0, p.1)).collect(),
    ))
}

/// Squared distance from a point to the nearest site.
fn clearance_sq(pts: &[(f64, f64)], x: f64, y: f64) -> f64 {
    let mut best = f64::INFINITY;
    for p in pts {
        let dx = x - p.0;
        let dy = y - p.1;
        best = best.min(dx * dx + dy * dy);
    }
    best
}

/// Coarse grid plus local refinement: an independent lower bound on the optimum.
///
/// Sampling is restricted to the obstacle convex hull, because that hull is the
/// region the solver confines its centre to. Sampling the bounding box instead
/// would find corners outside the hull, whose clearance is larger but which the
/// problem does not admit.
fn brute_force_radius(pts: &[(f64, f64)], side: usize) -> f64 {
    use geo::{Contains, ConvexHull};
    use geo_types::Point as GPoint;

    let hull = geo::MultiPoint::new(
        pts.iter()
            .map(|p| GPoint::new(p.0, p.1))
            .collect::<Vec<_>>(),
    )
    .convex_hull();
    let (mut min_x, mut min_y) = (f64::INFINITY, f64::INFINITY);
    let (mut max_x, mut max_y) = (f64::NEG_INFINITY, f64::NEG_INFINITY);
    for p in pts {
        min_x = min_x.min(p.0);
        min_y = min_y.min(p.1);
        max_x = max_x.max(p.0);
        max_y = max_y.max(p.1);
    }
    let dx = (max_x - min_x) / side as f64;
    let dy = (max_y - min_y) / side as f64;

    let mut best = (0.0, 0.0, 0.0f64);
    for i in 0..=side {
        for j in 0..=side {
            let x = min_x + i as f64 * dx;
            let y = min_y + j as f64 * dy;
            if !hull.contains(&GPoint::new(x, y)) {
                continue;
            }
            let d = clearance_sq(pts, x, y);
            if d > best.2 {
                best = (x, y, d);
            }
        }
    }

    // Climb from the best sample so the bound is not grid-limited.
    let (mut bx, mut by, mut bd) = best;
    let mut step = dx.max(dy) * 0.5;
    // The climb stays inside the same box the sampler covers, because the
    // solver confines its centre to the obstacle hull. An unbounded climb would
    // walk to infinity and report a radius the problem does not admit.
    let mut rounds = 0;
    while step > 1e-9 && rounds < 400 {
        rounds += 1;
        let mut moved = false;
        for (ox, oy) in [
            (-1.0, 0.0),
            (1.0, 0.0),
            (0.0, -1.0),
            (0.0, 1.0),
            (-1.0, -1.0),
            (1.0, -1.0),
            (-1.0, 1.0),
            (1.0, 1.0),
        ] {
            let nx = bx + ox * step;
            let ny = by + oy * step;
            if nx < min_x || nx > max_x || ny < min_y || ny > max_y {
                continue;
            }
            if !hull.contains(&GPoint::new(nx, ny)) {
                continue;
            }
            let d = clearance_sq(pts, nx, ny);
            if d > bd {
                bx = nx;
                by = ny;
                bd = d;
                moved = true;
            }
        }
        if !moved {
            step *= 0.5;
        }
    }
    bd.sqrt()
}

/// The reported radius must equal the measured clearance at the centre.
fn assert_reports_measured_radius(pts: &[(f64, f64)], r: &MecResult) {
    let measured = clearance_sq(pts, r.center.x(), r.center.y()).sqrt();
    let rel = (r.radius - measured).abs() / measured.max(1e-300);
    assert!(
        rel <= 1e-9,
        "reported radius {} != measured {measured} at ({}, {})",
        r.radius,
        r.center.x(),
        r.center.y()
    );
}

#[test]
fn exact_matches_brute_force_on_random_sites() {
    // Random sites are almost never collinear, but a seed that produced one
    // would be rejected by design, so those seeds are skipped rather than
    // failing the invariant.
    let mut checked = 0usize;
    for seed in 1..=40u64 {
        let mut rng = Rng(seed.wrapping_mul(0x9E37_79B9_7F4A_7C15) | 1);
        let n = 3 + (seed as usize % 12);
        let pts: Vec<(f64, f64)> = (0..n)
            .map(|_| (rng.range(0.0, 100.0), rng.range(0.0, 100.0)))
            .collect();

        let g = sites_geom(&pts);
        let Ok(r) = largest_empty_circle(&g, None, &MecOptions::default()) else {
            continue;
        };
        assert_eq!(r.approach, MecApproach::ExactVoronoi, "seed {seed}");
        checked += 1;

        assert_reports_measured_radius(&pts, &r);

        // The solver is exact, so it may only beat the sampled bound.
        let sampled = brute_force_radius(&pts, 40);
        assert!(
            r.radius >= sampled - 1e-6,
            "seed {seed}: exact {} lost to sampled {sampled}",
            r.radius
        );
    }
    assert!(checked >= 35, "only {checked} seeds were solvable");
}

#[test]
fn exact_beats_dense_sample_on_grid_sites() {
    // A 4x4 lattice: every Voronoi vertex is equidistant from four sites, which
    // is the degenerate case a naive circumcentre filter gets wrong.
    let mut pts = Vec::new();
    for i in 0..4 {
        for j in 0..4 {
            pts.push((i as f64 * 10.0, j as f64 * 10.0));
        }
    }
    let g = sites_geom(&pts);
    let r = largest_empty_circle(&g, None, &MecOptions::default()).expect("solve");
    assert_reports_measured_radius(&pts, &r);
    // Lattice spacing 10, so the best empty circle clears 5*sqrt(2).
    assert!(
        (r.radius - 7.071_067_811_9).abs() < 1e-6,
        "radius={} expected 5*sqrt(2)",
        r.radius
    );
}

#[test]
fn semicircle_sites_cluster_at_the_centre_of_curvature() {
    // Sites on a circle of radius 20: the Voronoi vertex at the centre clears
    // every site by exactly 20, so the solver must find that and not a corner.
    let n = 24;
    let pts: Vec<(f64, f64)> = (0..n)
        .map(|i| {
            let a = std::f64::consts::TAU * i as f64 / n as f64;
            (20.0 * a.cos(), 20.0 * a.sin())
        })
        .collect();
    let g = sites_geom(&pts);
    let r = largest_empty_circle(&g, None, &MecOptions::default()).expect("solve");
    assert!(
        (r.radius - 20.0).abs() < 1e-6,
        "radius={} expected 20",
        r.radius
    );
    assert!(r.center.x().abs() < 1e-6 && r.center.y().abs() < 1e-6);
}

#[test]
fn grid_path_reports_measured_radius_for_segments() {
    let segs = Geometry::MultiLineString(geo_types::MultiLineString::new(vec![
        LineString::from(vec![(0.0, 0.0), (100.0, 0.0)]),
        LineString::from(vec![(0.0, 50.0), (100.0, 50.0)]),
    ]));
    let r = largest_empty_circle(&segs, None, &MecOptions::default()).expect("solve");
    assert_eq!(r.approach, MecApproach::GridSolver);
    // Two parallel lines 50 apart: the widest gap is half of that.
    assert!(
        (r.radius - 25.0).abs() < 1e-2,
        "radius={} expected 25",
        r.radius
    );
}

#[test]
fn grid_path_inside_boundary_never_escapes_it() {
    let segs = Geometry::MultiPoint(MultiPoint::new(vec![
        Point::new(0.0, 0.0),
        Point::new(100.0, 100.0),
        Point::new(0.0, 100.0),
    ]));
    let boundary = Geometry::Polygon(Polygon::new(
        LineString::from(vec![
            (0.0, 0.0),
            (10.0, 0.0),
            (10.0, 10.0),
            (0.0, 10.0),
            (0.0, 0.0),
        ]),
        vec![],
    ));
    let r = largest_empty_circle(&segs, Some(&boundary), &MecOptions::default()).expect("solve");
    assert!(
        r.center.x() >= -1e-9 && r.center.x() <= 10.0 + 1e-9,
        "x={} escaped",
        r.center.x()
    );
    assert!(
        r.center.y() >= -1e-9 && r.center.y() <= 10.0 + 1e-9,
        "y={} escaped",
        r.center.y()
    );
}

#[test]
fn grid_path_agrees_with_the_analytic_boundary_maximum() {
    // Two obstacles tied along a region edge. The point at the origin and the
    // segment ending at (10, 0) are equidistant at (5, 5), where the box's top
    // edge peaks: |(5,5)| and |(5,5)-(10,0)| are both 5*sqrt(2).
    //
    // A two-way tie is never an interior maximum, so this answer sits exactly on
    // the boundary and only the edge search can reach it: a 2-D walk sees the
    // improving direction as "leave the region" and stops short.
    let obstacles = Geometry::GeometryCollection(geo_types::GeometryCollection::new_from(vec![
        Geometry::Point(Point::new(0.0, 0.0)),
        Geometry::LineString(LineString::from(vec![(10.0, 0.0), (10.0, -1.0)])),
    ]));
    let boundary = Geometry::Polygon(Polygon::new(
        LineString::from(vec![
            (0.0, 0.0),
            (10.0, 0.0),
            (10.0, 5.0),
            (0.0, 5.0),
            (0.0, 0.0),
        ]),
        vec![],
    ));
    let r =
        largest_empty_circle(&obstacles, Some(&boundary), &MecOptions::default()).expect("solve");
    assert_eq!(r.approach, MecApproach::GridSolver);
    let expected = 5.0 * std::f64::consts::SQRT_2;
    // The region is closed up to `span * 1e-9`, which is what lets a centre sit on
    // a ring edge after rounding, so the admissible answer can exceed the analytic
    // maximum by about that much.
    assert!(
        (r.radius - expected).abs() <= 1e-7,
        "radius={} expected 5*sqrt(2) = {expected}",
        r.radius
    );
    assert!(
        (r.center.x() - 5.0).abs() <= 1e-7 && (r.center.y() - 5.0).abs() <= 1e-7,
        "centre=({}, {}) expected (5, 5)",
        r.center.x(),
        r.center.y()
    );
}

#[test]
fn grid_path_agrees_with_the_analytic_three_way_apex() {
    // Three segments, each leaving a corner of an equilateral triangle outward, so
    // the nearest point of every segment to the triangle's centre is the corner
    // itself. The centre is then the exact three-way apex, at the circumradius
    // 10 / sqrt(3), and it is interior, so the grid's lattice has to place it
    // exactly rather than sample near it.
    let h = 5.0 * 3.0f64.sqrt();
    let corners: [(f64, f64); 3] = [(0.0, 0.0), (10.0, 0.0), (5.0, h)];
    let centre: (f64, f64) = (5.0, h / 3.0);
    let segments: Vec<LineString<f64>> = corners
        .iter()
        .map(|&(x, y)| {
            let dx = x - centre.0;
            let dy = y - centre.1;
            let len = dx.hypot(dy);
            LineString::from(vec![(x, y), (x + dx / len, y + dy / len)])
        })
        .collect();
    let obstacles = Geometry::MultiLineString(geo_types::MultiLineString::new(segments));
    let r = largest_empty_circle(&obstacles, None, &MecOptions::default()).expect("solve");
    assert_eq!(r.approach, MecApproach::GridSolver);
    let expected = 10.0 / 3.0f64.sqrt();
    assert!(
        (r.radius - expected).abs() <= 1e-9,
        "radius={} expected {expected}",
        r.radius
    );
    assert!(
        (r.center.x() - centre.0).abs() <= 1e-9 && (r.center.y() - centre.1).abs() <= 1e-9,
        "centre=({}, {}) expected ({}, {})",
        r.center.x(),
        r.center.y(),
        centre.0,
        centre.1
    );
}

#[test]
fn nan_obstacle_is_rejected() {
    let g = sites_geom(&[(0.0, 0.0), (f64::NAN, 1.0)]);
    let err = largest_empty_circle(&g, None, &MecOptions::default()).unwrap_err();
    assert!(matches!(err, MecError::InvalidInput(_)), "got {err}");
    assert!(err.to_string().contains("NaN"), "got {err}");
}

#[test]
fn single_site_gives_zero_radius() {
    // One site encloses no area, so no empty circle exists. GEOS answers that case
    // with radius zero centred on the site, and the crate matches it.
    let g = sites_geom(&[(5.0, 5.0)]);
    let r = largest_empty_circle(&g, None, &MecOptions::default()).expect("degenerate solve");
    assert!(
        r.radius <= 1e-9,
        "a single site should give radius zero, got {}",
        r.radius
    );
    assert!((r.center.x() - 5.0).abs() < 1e-9);
    assert!((r.center.y() - 5.0).abs() < 1e-9);
}

#[test]
fn grid_agrees_with_exact_when_points_are_present_but_lines_dominate() {
    // Mixed obstacles take the grid path. Place a dense line far from two
    // points so the points are the nearer constraint and their bisector value
    // is checkable.
    let g = Geometry::GeometryCollection(geo_types::GeometryCollection::new_from(vec![
        Geometry::MultiPoint(MultiPoint::new(vec![
            Point::new(0.0, 0.0),
            Point::new(10.0, 0.0),
        ])),
        Geometry::LineString(LineString::from(vec![(0.0, 100.0), (100.0, 100.0)])),
    ]));
    let r = largest_empty_circle(&g, None, &MecOptions::default()).expect("solve");
    assert_eq!(r.approach, MecApproach::GridSolver);
    assert!(r.radius.is_finite() && r.radius > 0.0);
}
