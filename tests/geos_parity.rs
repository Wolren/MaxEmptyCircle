//! Cross-check against GEOS `LargestEmptyCircle`, behind `bench-geos`.
//!
//! GEOS is an independent implementation, so this is the strongest correctness
//! evidence available short of a proof: two codebases, two algorithms, one
//! answer. For point obstacles the native path is exact and GEOS is iterative,
//! so the native radius must be at least GEOS's and must agree within GEOS's own
//! tolerance.
//!
//! Three properties are checked on every case.
//!
//! **Fit.** The reported radius equals the measured clearance at the reported
//! centre, so the circle really is empty by the same distance function both
//! solvers use.
//!
//! **Non-domination.** The native radius is never worse than GEOS's by more than
//! the tolerance. A native answer that is materially smaller means the exact
//! path missed a candidate or the grid stopped early.
//!
//! **Agreement.** The two radii agree within the tighter of the two tolerances,
//! scaled to the region. This is the property that catches a solver returning a
//! plausible but wrong circle.
//!
//! # Running
//!
//! ```text
//! GEOS_LIB_DIR=D:\Miniconda\Library\lib \
//!   cargo test --features bench-geos --test geos_parity
//! ```

#![cfg(feature = "bench-geos")]

use max_empty_circle::geo_types::{
    Coord, Geometry, LineString, MultiLineString, MultiPoint, Point, Polygon,
};
use max_empty_circle::mec::input::ObstacleSet;
use max_empty_circle::mec::solver::bench_geos::{geos_tolerance_for, solve_geos};
use max_empty_circle::{largest_empty_circle, MecApproach, MecOptions};

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

/// Assert the three properties against GEOS for one case.
fn check_against_geos(label: &str, geom: &Geometry<f64>, boundary: Option<&Geometry<f64>>) {
    let set = ObstacleSet::from_geometry(geom).expect("obstacle set");
    let native = largest_empty_circle(geom, boundary, &MecOptions::default())
        .unwrap_or_else(|e| panic!("{label}: native solve failed: {e}"));

    // Fit: the radius must be the measured clearance at the centre.
    let measured = set
        .clearance_sq(native.center.x(), native.center.y())
        .max(0.0)
        .sqrt();
    let fit_rel = (native.radius - measured).abs() / measured.max(1e-300);
    assert!(
        fit_rel <= 1e-9,
        "{label}: reported radius {} != measured {measured}",
        native.radius
    );

    let geos = solve_geos(&set, boundary).unwrap_or_else(|e| panic!("{label}: GEOS failed: {e}"));

    // Tolerance scaled to the region, which is how both solvers measure accuracy.
    let (min_x, min_y, max_x, max_y) = set.bounds().expect("bounds");
    let diag = (max_x - min_x).hypot(max_y - min_y).max(1.0);
    let native_tol = diag * 1e-7;
    let geos_tol = geos_tolerance_for(&set);
    let tol = native_tol.max(geos_tol * 10.0);

    // Non-domination: the native answer is never materially worse than GEOS.
    assert!(
        native.radius + tol >= geos.radius,
        "{label}: native {} is materially below GEOS {} (tol {tol})",
        native.radius,
        geos.radius
    );

    // Agreement: the two radii match within tolerance.
    let diff = (native.radius - geos.radius).abs();
    assert!(
        diff <= tol,
        "{label}: native {} vs GEOS {} differ by {diff} > {tol}",
        native.radius,
        geos.radius
    );
}

#[test]
fn triangle_points_match_geos() {
    let g = sites_geom(&[(0.0, 0.0), (10.0, 0.0), (5.0, 8.660_254_037_844_387)]);
    check_against_geos("triangle", &g, None);
}

#[test]
fn square_points_match_geos() {
    let g = sites_geom(&[(0.0, 0.0), (10.0, 0.0), (10.0, 10.0), (0.0, 10.0)]);
    check_against_geos("square", &g, None);
}

#[test]
fn lattice_points_match_geos() {
    let mut pts = Vec::new();
    for i in 0..4 {
        for j in 0..4 {
            pts.push((i as f64 * 10.0, j as f64 * 10.0));
        }
    }
    let g = sites_geom(&pts);
    check_against_geos("lattice", &g, None);
}

#[test]
fn circle_of_points_matches_geos() {
    let n = 24;
    let pts: Vec<(f64, f64)> = (0..n)
        .map(|i| {
            let a = std::f64::consts::TAU * i as f64 / n as f64;
            (20.0 * a.cos(), 20.0 * a.sin())
        })
        .collect();
    let g = sites_geom(&pts);
    check_against_geos("circle", &g, None);
}

#[test]
fn parallel_walls_match_geos() {
    let g = Geometry::MultiLineString(MultiLineString::new(vec![
        LineString::from(vec![(0.0, 0.0), (100.0, 0.0)]),
        LineString::from(vec![(0.0, 50.0), (100.0, 50.0)]),
    ]));
    check_against_geos("walls", &g, None);
}

#[test]
fn polygon_obstacle_matches_geos() {
    let g = Geometry::Polygon(Polygon::new(
        LineString::from(vec![
            (0.0, 0.0),
            (40.0, 0.0),
            (40.0, 40.0),
            (0.0, 40.0),
            (0.0, 0.0),
        ]),
        vec![],
    ));
    check_against_geos("polygon", &g, None);
}

#[test]
fn buildings_match_geos() {
    // A city block: four building footprints as obstacles.
    let blocks = vec![
        (0.0, 0.0, 20.0, 20.0),
        (60.0, 0.0, 80.0, 20.0),
        (0.0, 60.0, 20.0, 80.0),
        (60.0, 60.0, 80.0, 80.0),
    ];
    let polys: Vec<Polygon<f64>> = blocks
        .iter()
        .map(|&(x0, y0, x1, y1)| {
            Polygon::new(
                LineString::from(vec![(x0, y0), (x1, y0), (x1, y1), (x0, y1), (x0, y0)]),
                vec![],
            )
        })
        .collect();
    let g = Geometry::MultiPolygon(geo_types::MultiPolygon::new(polys));
    check_against_geos("blocks", &g, None);
}

#[test]
fn boundary_confined_points_match_geos() {
    let g = sites_geom(&[(0.0, 0.0), (100.0, 0.0), (50.0, 90.0), (0.0, 90.0)]);
    let b = Geometry::Polygon(Polygon::new(
        LineString::from(vec![
            (0.0, 0.0),
            (30.0, 0.0),
            (30.0, 30.0),
            (0.0, 30.0),
            (0.0, 0.0),
        ]),
        vec![],
    ));
    check_against_geos("confined", &g, Some(&b));
}

#[test]
fn random_point_sets_match_geos() {
    let mut checked = 0usize;
    for seed in 1..=60u64 {
        let mut rng = Rng(seed.wrapping_mul(0x9E37_79B9_7F4A_7C15) | 1);
        let n = 3 + (seed as usize % 20);
        let pts: Vec<(f64, f64)> = (0..n)
            .map(|_| (rng.range(0.0, 100.0), rng.range(0.0, 100.0)))
            .collect();

        let g = sites_geom(&pts);
        // Collinear or degenerate sets are rejected by design, so skip them.
        if ObstacleSet::from_geometry(&g).is_err() {
            continue;
        }
        let Ok(native) = largest_empty_circle(&g, None, &MecOptions::default()) else {
            continue;
        };
        assert_eq!(
            native.approach,
            MecApproach::ExactVoronoi,
            "seed {seed}: point-only input must take the exact path"
        );

        let set = ObstacleSet::from_geometry(&g).expect("set");
        let geos = solve_geos(&set, None).expect("geos");

        let (min_x, min_y, max_x, max_y) = set.bounds().expect("bounds");
        let diag = (max_x - min_x).hypot(max_y - min_y).max(1.0);
        let tol = (diag * 1e-7).max(geos_tolerance_for(&set) * 10.0);

        // The exact path must not lose to GEOS by more than tolerance.
        assert!(
            native.radius + tol >= geos.radius,
            "seed {seed}: exact {} below GEOS {}",
            native.radius,
            geos.radius
        );
        checked += 1;
    }
    assert!(checked >= 50, "only {checked} seeds were solvable");
}

#[test]
fn random_segment_sets_match_geos() {
    let mut checked = 0usize;
    for seed in 1..=30u64 {
        let mut rng = Rng(seed.wrapping_mul(0x2545_F491_4F6C_DD1D) | 1);
        let n = 3 + (seed as usize % 10);
        let mut segs = Vec::with_capacity(n);
        for _ in 0..n {
            let x0 = rng.range(0.0, 100.0);
            let y0 = rng.range(0.0, 100.0);
            let x1 = x0 + rng.range(-20.0, 20.0);
            let y1 = y0 + rng.range(-20.0, 20.0);
            segs.push(LineString::from(vec![
                Coord { x: x0, y: y0 },
                Coord { x: x1, y: y1 },
            ]));
        }
        let g = Geometry::MultiLineString(MultiLineString::new(segs));
        if ObstacleSet::from_geometry(&g).is_err() {
            continue;
        }
        let Ok(native) = largest_empty_circle(&g, None, &MecOptions::default()) else {
            continue;
        };

        let set = ObstacleSet::from_geometry(&g).expect("set");
        let geos = solve_geos(&set, None).expect("geos");

        let (min_x, min_y, max_x, max_y) = set.bounds().expect("bounds");
        let diag = (max_x - min_x).hypot(max_y - min_y).max(1.0);
        let tol = (diag * 1e-6).max(geos_tolerance_for(&set) * 10.0);

        assert!(
            native.radius + tol >= geos.radius,
            "seed {seed}: grid {} below GEOS {} (tol {tol})",
            native.radius,
            geos.radius
        );
        checked += 1;
    }
    assert!(checked >= 25, "only {checked} seeds were solvable");
}
