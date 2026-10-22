//! Complexity regression guard: a superlinear path must not reappear unnoticed.
//!
//! The solve time is measured at doubling input sizes and the growth exponent is
//! fitted on the log-log points. A linear solver fits near 1.0, a quadratic one
//! near 2.0, a cubic one near 3.0.
//!
//! The bound below is deliberately loose. A timing assertion tight enough to catch
//! a 20 percent wobble would fail on a loaded CI machine and get muted, while an
//! accidental quadratic path moves the exponent far past it. The check exists to
//! catch a complexity class change, not a constant-factor drift.
#![cfg(not(debug_assertions))]

use std::time::Instant;

use max_empty_circle::geo_types::{
    Coord, Geometry, LineString, MultiLineString, MultiPoint, Point,
};
use max_empty_circle::{largest_empty_circle, MecOptions};

/// Deterministic xorshift, so a failure reproduces from the seed alone.
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

/// Time one solve class across sizes, returning `(n, milliseconds)` pairs.
fn sweep(sizes: &[usize], build: impl Fn(usize) -> Geometry<f64>) -> Vec<(usize, f64)> {
    let opts = MecOptions::default();
    let mut out = Vec::with_capacity(sizes.len());
    for &n in sizes {
        let geom = build(n);
        let start = Instant::now();
        let r = largest_empty_circle(&geom, None, &opts).expect("solve");
        let ms = start.elapsed().as_secs_f64() * 1e3;
        std::hint::black_box(r.radius);
        out.push((n, ms));
    }
    out
}

/// Exponent of the fitted `t = c * n^k`.
fn exponent(points: &[(usize, f64)]) -> f64 {
    let usable: Vec<(f64, f64)> = points
        .iter()
        .filter(|(n, ms)| *n > 0 && *ms > 0.0)
        .map(|(n, ms)| ((*n as f64).ln(), ms.ln()))
        .collect();
    if usable.len() < 2 {
        return 0.0;
    }
    let count = usable.len() as f64;
    let sum_x: f64 = usable.iter().map(|p| p.0).sum();
    let sum_y: f64 = usable.iter().map(|p| p.1).sum();
    let sum_xx: f64 = usable.iter().map(|p| p.0 * p.0).sum();
    let sum_xy: f64 = usable.iter().map(|p| p.0 * p.1).sum();
    let denom = count * sum_xx - sum_x * sum_x;
    if denom.abs() < f64::EPSILON {
        return 0.0;
    }
    (count * sum_xy - sum_x * sum_y) / denom
}

fn points_geom(n: usize) -> Geometry<f64> {
    let mut rng = Rng(0x9E37_79B9_7F4A_7C15 ^ n as u64);
    Geometry::MultiPoint(MultiPoint::new(
        (0..n)
            .map(|_| Point::new(rng.range(0.0, 1000.0), rng.range(0.0, 1000.0)))
            .collect::<Vec<_>>(),
    ))
}

fn segments_geom(n: usize) -> Geometry<f64> {
    let mut rng = Rng(0x2545_F491_4F6C_DD1D ^ n as u64);
    Geometry::MultiLineString(MultiLineString::new(
        (0..n)
            .map(|_| {
                let x = rng.range(0.0, 1000.0);
                let y = rng.range(0.0, 1000.0);
                LineString::from(vec![
                    Coord { x, y },
                    Coord {
                        x: x + rng.range(-40.0, 40.0),
                        y: y + rng.range(-40.0, 40.0),
                    },
                ])
            })
            .collect::<Vec<_>>(),
    ))
}

#[test]
fn point_growth_stays_subquadratic() {
    // Above the exact budget the grid answers, so this measures the grid path.
    let points = sweep(&[400, 1_600, 6_400, 25_600], points_geom);
    let k = exponent(&points);
    assert!(
        k < 1.8,
        "point solve grows at n^{k:.2}, which is past the subquadratic bound: {points:?}"
    );
}

#[test]
fn segment_growth_stays_subquadratic() {
    let points = sweep(&[250, 1_000, 4_000], segments_geom);
    let k = exponent(&points);
    assert!(
        k < 1.8,
        "segment solve grows at n^{k:.2}, which is past the subquadratic bound: {points:?}"
    );
}

#[test]
fn exact_path_stays_subquadratic_below_budget() {
    // The exact path is Delaunay plus linear boundary families, so it should be
    // close to n log n rather than the n^3 the triple loop used to be.
    let points = sweep(&[25, 50, 100], points_geom);
    let k = exponent(&points);
    assert!(
        k < 2.5,
        "exact path grows at n^{k:.2}, which is past the n log n expectation: {points:?}"
    );
}

#[test]
fn large_input_completes() {
    // A guard that the solver does not degrade into an unbounded loop at scale.
    let start = Instant::now();
    let geom = points_geom(100_000);
    let r = largest_empty_circle(&geom, None, &MecOptions::default()).expect("solve");
    let elapsed = start.elapsed().as_secs_f64();
    assert!(r.radius > 0.0);
    assert!(
        elapsed < 120.0,
        "100k points took {elapsed:.1}s, which suggests a complexity regression"
    );
}
