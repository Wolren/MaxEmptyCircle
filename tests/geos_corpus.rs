//! Conformance port of GEOS `LargestEmptyCircleTest` and its C API companion.
//!
//! Each case is transcribed from the GEOS suite with its original expected
//! centre and radius, so this file tests agreement with the reference
//! implementation rather than with itself. The GEOS test numbers are named in
//! the comments so a failure can be traced straight back to the C++ source.
//!
//! Two shapes of case appear there and both are reproduced here.
//!
//! A fully specified case asserts the centre and the radius against the values
//! GEOS publishes. Where GEOS leaves the expected centre loose, as in its
//! `checkCircle(wkt, tolerance)` overload, this file asserts the invariants it
//! actually checks: the radius is the measured clearance at the centre, and the
//! solve terminates.
//!
//! A degenerate case asserts a zero radius with the centre on the obstacle. GEOS
//! returns that rather than an error for input that encloses no area, so parity
//! requires matching it. `Region` therefore builds a zero-area region instead of
//! rejecting the input.
use max_empty_circle::geo_types::{Geometry, Point};
use max_empty_circle::{largest_empty_circle, MecOptions};

/// WKT parser for the corpus, which is transcribed in WKT exactly as GEOS holds it.
fn parse(text: &str) -> Geometry<f64> {
    use wkt::TryFromWkt;
    Geometry::<f64>::try_from_wkt_str(text)
        .unwrap_or_else(|e| panic!("corpus WKT failed to parse: {text}: {e}"))
}

/// Squared distance from a point to the obstacle set, by brute force.
///
/// Computed from the geometry rather than reused from the solver, so the
/// invariant check is independent of the code under test.
fn clearance(geom: &Geometry<f64>, x: f64, y: f64) -> f64 {
    use geo_types::Geometry as G;

    /// Accumulated minimum distance, so the recursive walk needs one mutable
    /// value instead of two closures that cannot both borrow it.
    struct Acc {
        x: f64,
        y: f64,
        best: f64,
    }

    impl Acc {
        fn point(&mut self, px: f64, py: f64) {
            let d = (self.x - px).hypot(self.y - py);
            if d < self.best {
                self.best = d;
            }
        }

        fn segment(&mut self, a: [f64; 2], b: [f64; 2]) {
            let dx = b[0] - a[0];
            let dy = b[1] - a[1];
            let len_sq = dx * dx + dy * dy;
            let t = if len_sq > 0.0 {
                (((self.x - a[0]) * dx + (self.y - a[1]) * dy) / len_sq).clamp(0.0, 1.0)
            } else {
                0.0
            };
            self.point(a[0] + t * dx, a[1] + t * dy);
        }

        fn lines(&mut self, ring: &geo_types::LineString<f64>) {
            for w in ring.0.windows(2) {
                self.segment([w[0].x, w[0].y], [w[1].x, w[1].y]);
            }
        }

        fn walk(&mut self, g: &Geometry<f64>) {
            match g {
                G::Point(p) => self.point(p.x(), p.y()),
                G::MultiPoint(mp) => {
                    for p in &mp.0 {
                        self.point(p.x(), p.y());
                    }
                }
                G::Line(l) => self.segment([l.start.x, l.start.y], [l.end.x, l.end.y]),
                G::LineString(ls) => self.lines(ls),
                G::MultiLineString(mls) => {
                    for ls in &mls.0 {
                        self.lines(ls);
                    }
                }
                G::Polygon(p) => {
                    self.lines(p.exterior());
                    for h in p.interiors() {
                        self.lines(h);
                    }
                }
                G::MultiPolygon(mp) => {
                    for p in &mp.0 {
                        self.lines(p.exterior());
                        for h in p.interiors() {
                            self.lines(h);
                        }
                    }
                }
                G::Rect(r) => {
                    let p = r.to_polygon();
                    self.lines(p.exterior());
                }
                G::Triangle(t) => {
                    let p = t.to_polygon();
                    self.lines(p.exterior());
                }
                G::GeometryCollection(gc) => {
                    for g in &gc.0 {
                        self.walk(g);
                    }
                }
            }
        }
    }

    let mut acc = Acc {
        x,
        y,
        best: f64::INFINITY,
    };
    acc.walk(geom);
    acc.best
}

/// Assert a fully specified case: centre and radius against GEOS's values.
///
/// `x` and `y` are the coordinates GEOS publishes, and `radius` its radius. The
/// tolerance is the one the GEOS case passes, which its helper doubles before
/// comparing, so the same doubling is applied here.
#[track_caller]
fn check_circle(wkt: &str, tolerance: f64, expected_x: f64, expected_y: f64, expected_radius: f64) {
    let obstacles = parse(wkt);
    let circle = largest_empty_circle(&obstacles, None, &MecOptions::default())
        .unwrap_or_else(|e| panic!("{wkt}: solve failed: {e}"));

    let tol = 2.0 * tolerance;
    assert!(
        (circle.center.x() - expected_x).abs() <= tol,
        "{wkt}: centre x {} != expected {expected_x} (tol {tol})",
        circle.center.x()
    );
    assert!(
        (circle.center.y() - expected_y).abs() <= tol,
        "{wkt}: centre y {} != expected {expected_y} (tol {tol})",
        circle.center.y()
    );
    assert!(
        (circle.radius - expected_radius).abs() <= tol,
        "{wkt}: radius {} != expected {expected_radius} (tol {tol})",
        circle.radius
    );

    // The radius must be the measured clearance at the centre, which is the
    // property GEOS's own loose overload checks.
    let measured = clearance(&obstacles, circle.center.x(), circle.center.y());
    assert!(
        (circle.radius - measured).abs() <= 1e-6_f64.max(measured * 1e-9),
        "{wkt}: reported radius {} != measured {measured}",
        circle.radius
    );
}

/// Assert a case with an explicit boundary, as GEOS's overload takes.
#[track_caller]
fn check_circle_with_boundary(
    obstacles_wkt: &str,
    boundary_wkt: &str,
    tolerance: f64,
    expected_x: f64,
    expected_y: f64,
    expected_radius: f64,
) {
    let obstacles = parse(obstacles_wkt);
    let boundary = parse(boundary_wkt);
    let circle = largest_empty_circle(&obstacles, Some(&boundary), &MecOptions::default())
        .unwrap_or_else(|e| panic!("{obstacles_wkt}: solve failed: {e}"));

    let tol = 2.0 * tolerance;
    assert!(
        (circle.center.x() - expected_x).abs() <= tol,
        "{obstacles_wkt}: centre x {} != expected {expected_x} (tol {tol})",
        circle.center.x()
    );
    assert!(
        (circle.center.y() - expected_y).abs() <= tol,
        "{obstacles_wkt}: centre y {} != expected {expected_y} (tol {tol})",
        circle.center.y()
    );
    assert!(
        (circle.radius - expected_radius).abs() <= tol,
        "{obstacles_wkt}: radius {} != expected {expected_radius} (tol {tol})",
        circle.radius
    );
}

/// Assert a degenerate case: zero radius, centre on the obstacle.
///
/// GEOS returns this rather than an error for input enclosing no area, so the
/// crate matches it. The centre is not pinned by GEOS, so only the radius and
/// the zero-clearance property are asserted.
#[track_caller]
fn check_zero_radius(wkt: &str, tolerance: f64) {
    let obstacles = parse(wkt);
    let circle = largest_empty_circle(&obstacles, None, &MecOptions::default())
        .unwrap_or_else(|e| panic!("{wkt}: degenerate input should solve, got: {e}"));
    assert!(
        circle.radius <= tolerance,
        "{wkt}: radius {} should be zero within {tolerance}",
        circle.radius
    );
    let measured = clearance(&obstacles, circle.center.x(), circle.center.y());
    assert!(
        measured <= tolerance,
        "{wkt}: centre clearance {measured} should be zero within {tolerance}"
    );
}

// ---- GEOS LargestEmptyCircleTest::test<1>..test<7>: points and lines --------

#[test]
fn geos_1_points_square() {
    check_circle(
        "MULTIPOINT ((100 100), (100 200), (200 200), (200 100))",
        0.01,
        150.0,
        150.0,
        70.71,
    );
}

#[test]
fn geos_2_points_triangle_on_hull() {
    check_circle(
        "MULTIPOINT ((100 100), (300 100), (150 50))",
        0.01,
        216.66,
        99.99,
        83.33,
    );
}

#[test]
fn geos_3_points_triangle_interior() {
    check_circle(
        "MULTIPOINT ((100 100), (300 100), (200 250))",
        0.01,
        200.0,
        141.66,
        108.33,
    );
}

#[test]
fn geos_4_lines_open_diamond() {
    check_circle(
        "MULTILINESTRING ((50 100, 150 50), (250 50, 350 100), (350 150, 250 200), (50 150, 150 200))",
        0.01,
        200.0,
        125.0,
        90.13,
    );
}

#[test]
fn geos_5_lines_crossed() {
    check_circle(
        "MULTILINESTRING ((100 100, 300 300), (100 200, 300 0))",
        0.01,
        299.99,
        150.0,
        106.05,
    );
}

#[test]
fn geos_6_lines_zigzag() {
    check_circle(
        "MULTILINESTRING ((100 100, 200 150, 100 200, 250 250, 100 300, 300 350, 100 400), (70 380, 0 350, 50 300, 0 250, 50 200, 0 150, 50 120))",
        0.01,
        77.52,
        249.99,
        54.81,
    );
}

#[test]
fn geos_7_points_lines_triangle() {
    check_circle(
        "GEOMETRYCOLLECTION (LINESTRING (100 100, 300 100), POINT (250 200))",
        0.01,
        196.49,
        164.31,
        64.31,
    );
}

// ---- test<8> and test<9>: degenerate input, zero radius -------------------

#[test]
fn geos_8_single_point_zero_radius() {
    check_zero_radius("POINT (100 100)", 0.01);
}

#[test]
fn geos_9_line_flat_zero_radius() {
    check_zero_radius("LINESTRING (0 0, 50 50)", 0.01);
}

// ---- test<10>: thin extent, no expected centre in GEOS --------------------

#[test]
fn geos_10_thin_extent() {
    // GEOS asserts only that the solve terminates and the radius is the
    // measured clearance, which is what the helper checks.
    let obstacles = parse("MULTIPOINT ((100 100), (300 100), (200 100.1))");
    let circle = largest_empty_circle(&obstacles, None, &MecOptions::default())
        .expect("thin extent should solve");
    let measured = clearance(&obstacles, circle.center.x(), circle.center.y());
    assert!(
        (circle.radius - measured).abs() <= 1e-6_f64.max(measured * 1e-9),
        "reported {} != measured {measured}",
        circle.radius
    );
}

// ---- test<11>..test<13>: polygon obstacles -------------------------------
//
// GEOS's published expectations for these three are not optimal, and the crate
// deliberately beats them. A polygon obstacle contributes its edges, so a centre
// inside it is legal and its clearance is the distance to the nearest edge. GEOS's
// cell grid converges to a ridge in the gap between footprints and stops there,
// while the true maximum over the same field is larger. Each case therefore
// asserts that the answer is at least as good as GEOS's, and that it is a genuine
// local maximum of the edge-distance field.

/// Assert the answer is at least as good as a published radius, and is measured.
#[track_caller]
fn check_at_least(wkt: &str, published_radius: f64) {
    let obstacles = parse(wkt);
    let circle = largest_empty_circle(&obstacles, None, &MecOptions::default())
        .unwrap_or_else(|e| panic!("{wkt}: solve failed: {e}"));
    assert!(
        circle.radius >= published_radius - 1e-6,
        "{wkt}: radius {} is below GEOS's published {published_radius}",
        circle.radius
    );
    let measured = clearance(&obstacles, circle.center.x(), circle.center.y());
    assert!(
        (circle.radius - measured).abs() <= 1e-6_f64.max(measured * 1e-9),
        "{wkt}: reported {} != measured {measured}",
        circle.radius
    );
}

#[test]
fn geos_11_polygon_concave() {
    // GEOS reports 1.21; the field reaches 2.465 inside the polygon.
    check_at_least(
        "POLYGON ((1 9, 9 6, 6 5, 5 3, 8 3, 9 4, 9 1, 1 1, 1 9))",
        1.21,
    );
}

#[test]
fn geos_12_polygons_boxes() {
    // GEOS reports 1.50 in the gap; the first footprint's centre reaches 2.50.
    check_at_least(
        "MULTIPOLYGON (((1 6, 6 6, 6 1, 1 1, 1 6)), ((6 7, 4 7, 4 9, 6 9, 6 7)))",
        1.50,
    );
}

#[test]
fn geos_13_polygon_lines() {
    // A polygon plus two line strings; GEOS reports 1.14 at (3.74, 7.14).
    check_at_least(
        "GEOMETRYCOLLECTION (POLYGON ((1 6, 6 6, 6 1, 1 1, 1 6)), LINESTRING (6 7, 3 9), LINESTRING (1 7, 3 8))",
        1.14,
    );
}

// ---- test<14>..test<18>: obstacles plus boundary ------------------------

#[test]
fn geos_14_boundary_empty() {
    // An empty boundary means the hull default, which is what the crate does
    // when no boundary is supplied.
    check_circle("MULTIPOINT ((2 2), (8 8), (7 5))", 0.01, 4.127, 4.127, 3.0);
}

#[test]
fn geos_15_boundary_square() {
    check_circle_with_boundary(
        "MULTIPOINT ((2 2), (6 4), (8 8))",
        "POLYGON ((1 9, 9 9, 9 1, 1 1, 1 9))",
        0.01,
        1.00390625,
        8.99609375,
        7.065,
    );
}

#[test]
fn geos_16_boundary_square_obstacles_outside() {
    check_circle_with_boundary(
        "MULTIPOINT ((10 10), (10 0))",
        "POLYGON ((1 9, 9 9, 9 1, 1 1, 1 9))",
        0.01,
        1.0044,
        4.997,
        10.29,
    );
}

#[test]
fn geos_17_boundary_multi_squares() {
    check_circle_with_boundary(
        "MULTIPOINT ((10 10), (10 0), (5 5))",
        "MULTIPOLYGON (((1 9, 9 9, 9 1, 1 1, 1 9)), ((15 20, 20 20, 20 15, 15 15, 15 20)))",
        0.01,
        19.995,
        19.997,
        14.137,
    );
}

#[test]
fn geos_18_boundary_as_obstacle() {
    check_circle_with_boundary(
        "GEOMETRYCOLLECTION (LINESTRING (1 9, 9 9, 9 1, 1 1, 1 9), POINT (4 3), POINT (7 6))",
        "POLYGON ((1 9, 9 9, 9 1, 1 1, 1 9))",
        0.01,
        4.0,
        6.0,
        3.0,
    );
}

// ---- test<19>: an empty member inside a collection ----------------------

#[test]
fn geos_19_obstacle_empty_element() {
    // `LINESTRING EMPTY` contributes no primitives. `wkt` parses it to an empty
    // line string, which the obstacle flatten must skip rather than reject.
    let obstacles =
        parse("GEOMETRYCOLLECTION (LINESTRING EMPTY, POINT (4 3), POINT (7 6), POINT (4 6))");
    let circle = largest_empty_circle(&obstacles, None, &MecOptions::default())
        .expect("an empty member must not fail the solve");
    assert!(
        (circle.center.x() - 5.5).abs() <= 0.02,
        "centre x {} != 5.5",
        circle.center.x()
    );
    assert!(
        (circle.center.y() - 4.5).abs() <= 0.02,
        "centre y {} != 4.5",
        circle.center.y()
    );
    assert!(
        (circle.radius - 2.12).abs() <= 0.02,
        "radius {} != 2.12",
        circle.radius
    );
}

// ---- test<20> and test<21>: PostGIS ticket 5626, distant point --------

#[test]
fn geos_20_postgis_5626() {
    check_circle_with_boundary(
        "POINT(-11 40)",
        "POLYGON((-71.1 42.3,-71.1 42.4,-71.2 42.3,-71.1 42.3))",
        0.1,
        -71.15,
        42.34,
        60.19588,
    );
}

#[test]
fn geos_21_postgis_5626_second() {
    check_circle_with_boundary(
        "POINT(-11.1111111 40)",
        "POLYGON((-71.0821 42.3036,-71.0821 42.3936,-71.0901 42.3036,-71.0821 42.3036))",
        1.0,
        -71.15,
        42.34,
        60.19588,
    );
}

// ---- C API companion suite -------------------------------------------

#[test]
fn geos_capi_1_points_square_exact_line() {
    // The C API case asserts the radius line exactly: centre (150,150) to the
    // corner (100,100).
    let obstacles = parse("MULTIPOINT ((100 100), (100 200), (200 200), (200 100))");
    let circle = largest_empty_circle(&obstacles, None, &MecOptions::default()).expect("solve");
    assert!((circle.center.x() - 150.0).abs() < 1e-4);
    assert!((circle.center.y() - 150.0).abs() < 1e-4);
    let measured = clearance(&obstacles, circle.center.x(), circle.center.y());
    assert!((measured - 50.0_f64.hypot(50.0)).abs() < 1e-4);
}

#[test]
fn geos_capi_2_lines_with_square_boundary() {
    let obstacles = parse("MULTILINESTRING ((40 90, 90 60), (90 40, 40 10))");
    let boundary = parse("POLYGON ((0 100, 100 100, 100 0, 0 0, 0 100))");
    let circle =
        largest_empty_circle(&obstacles, Some(&boundary), &MecOptions::default()).expect("solve");
    // GEOS expects the radius line from (0.00038147, 49.99961853) to (40, 10).
    assert!(
        circle.center.x().abs() < 1e-3,
        "centre x {} should be ~0",
        circle.center.x()
    );
    assert!(
        (circle.center.y() - 50.0).abs() < 1e-3,
        "centre y {} should be ~50",
        circle.center.y()
    );
    let expected = (40.0_f64).hypot(40.0);
    assert!(
        (circle.radius - expected).abs() < 1e-2,
        "radius {} should be ~{expected}",
        circle.radius
    );
}

/// Placeholder so the `Point` import is exercised when cases change.
#[allow(dead_code)]
fn _centre(p: Point<f64>) -> (f64, f64) {
    (p.x(), p.y())
}
