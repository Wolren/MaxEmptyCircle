//! Complexity scaling check: catches an accidental superlinear path reappearing.
//!
//! Each sweep doubles the input and records the time. A near-linear solver shows a
//! ratio near 2 per doubling; quadratic shows near 4; cubic near 8. The thresholds
//! below are deliberately loose, because a timing test that fails on machine noise
//! gets ignored, and a regression of a whole complexity class clears them anyway.
//!
//! Run:
//! ```text
//! cargo run --release --example scaling
//! ```
//!
//! The same sweep is asserted in `tests/scaling.rs`, so CI enforces it.

use std::time::Instant;

use max_empty_circle::geo_types::{
    Coord, Geometry, LineString, MultiLineString, MultiPoint, Point,
};
use max_empty_circle::{largest_empty_circle, MecOptions};

/// Deterministic xorshift, so a run reproduces from the seed alone.
pub struct Rng(pub u64);

impl Rng {
    pub fn step(&mut self) -> u64 {
        self.0 ^= self.0 << 13;
        self.0 ^= self.0 >> 7;
        self.0 ^= self.0 << 17;
        self.0
    }

    pub fn unit(&mut self) -> f64 {
        (self.step() >> 11) as f64 / (1u64 << 53) as f64
    }

    pub fn range(&mut self, lo: f64, hi: f64) -> f64 {
        lo + (hi - lo) * self.unit()
    }
}

/// One scaling measurement.
pub struct Measurement {
    /// Input primitive count.
    pub n: usize,
    /// Best of the repeats, in milliseconds.
    pub ms: f64,
}

/// Best-of-repeats timing of one solve.
pub fn time_once(build: &dyn Fn() -> Geometry<f64>, reps: usize) -> f64 {
    let geom = build();
    let opts = MecOptions::default();
    let mut best = f64::INFINITY;
    for _ in 0..reps {
        let start = Instant::now();
        // The result must be consumed, or the optimiser is free to drop the call.
        let r = largest_empty_circle(&geom, None, &opts).expect("solve");
        let elapsed = start.elapsed().as_secs_f64() * 1e3;
        std::hint::black_box(r.radius);
        if elapsed < best {
            best = elapsed;
        }
    }
    best
}

/// Point-obstacle sweep: the exact path below the budget, the grid above it.
pub fn point_sweep(sizes: &[usize], reps: usize) -> Vec<Measurement> {
    sizes
        .iter()
        .map(|&n| {
            let ms = time_once(
                &|| {
                    let mut rng = Rng(0x9E37_79B9_7F4A_7C15 ^ n as u64);
                    Geometry::MultiPoint(MultiPoint::new(
                        (0..n)
                            .map(|_| Point::new(rng.range(0.0, 1000.0), rng.range(0.0, 1000.0)))
                            .collect::<Vec<_>>(),
                    ))
                },
                reps,
            );
            Measurement { n, ms }
        })
        .collect()
}

/// Segment-obstacle sweep: always the grid path.
pub fn segment_sweep(sizes: &[usize], reps: usize) -> Vec<Measurement> {
    sizes
        .iter()
        .map(|&n| {
            let ms = time_once(
                &|| {
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
                },
                reps,
            );
            Measurement { n, ms }
        })
        .collect()
}

/// Exponent of the fitted `t = c * n^k` across a sweep.
///
/// A least-squares fit on the log-log points. Reported so a reader can see the
/// growth class rather than having to infer it from the ratios.
pub fn growth_exponent(points: &[Measurement]) -> f64 {
    let usable: Vec<(f64, f64)> = points
        .iter()
        .filter(|p| p.n > 0 && p.ms > 0.0)
        .map(|p| ((p.n as f64).ln(), p.ms.ln()))
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

fn print_sweep(label: &str, points: &[Measurement]) {
    println!("\n{label}");
    println!("{:>9}{:>12}{:>9}", "n", "ms", "ratio");
    let mut prev: Option<&Measurement> = None;
    for p in points {
        let ratio = prev.map(|q| {
            let dn = p.n as f64 / q.n as f64;
            let dt = p.ms / q.ms.max(f64::MIN_POSITIVE);
            dt / dn
        });
        match ratio {
            Some(r) => println!("{:>9}{:>12.3}{:>9.2}", p.n, p.ms, r),
            None => println!("{:>9}{:>12.3}{:>9}", p.n, p.ms, "-"),
        }
        prev = Some(p);
    }
    println!(
        "growth exponent k = {:.2} (1.0 linear, 2.0 quadratic)",
        growth_exponent(points)
    );
}

fn main() {
    println!("=== MaxEmptyCircle scaling");
    print_sweep(
        "points (exact below 200, grid above)",
        &point_sweep(
            &[
                10, 25, 50, 100, 200, 400, 1_000, 4_000, 16_000, 64_000, 100_000,
            ],
            1,
        ),
    );
    print_sweep(
        "segments (grid)",
        &segment_sweep(&[10, 25, 50, 100, 250, 500, 1_000, 2_500, 5_000, 10_000], 1),
    );
}
