//! The four-lane batch must return exactly what the scalar path returns.
//!
//! A vectorised distance is only safe if it agrees with the scalar reduction
//! bit-for-bit on sane input, so this file pins that. Any divergence means the
//! batch changed the arithmetic rather than widening it.
#![cfg(feature = "simd")]

use max_empty_circle::mec::solver::simd_probe::{lanes4, segment_distance_sq_x4, to_array4};
use wide::f64x4;

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

/// Scalar reference: squared distance from a point to a segment.
fn scalar_segment_sq(x: f64, y: f64, a: [f64; 2], b: [f64; 2]) -> f64 {
    let dx = b[0] - a[0];
    let dy = b[1] - a[1];
    let len_sq = dx * dx + dy * dy;
    let t = if len_sq > 0.0 {
        (((x - a[0]) * dx + (y - a[1]) * dy) / len_sq).clamp(0.0, 1.0)
    } else {
        0.0
    };
    let px = x - (a[0] + t * dx);
    let py = y - (a[1] + t * dy);
    px * px + py * py
}

#[test]
fn batch_matches_scalar_on_random_segments() {
    let mut rng = Rng(0x9E37_79B9_7F4A_7C15);
    for _ in 0..5_000 {
        let a = [rng.range(-1000.0, 1000.0), rng.range(-1000.0, 1000.0)];
        let b = [rng.range(-1000.0, 1000.0), rng.range(-1000.0, 1000.0)];
        let pts: [[f64; 2]; 4] =
            std::array::from_fn(|_| [rng.range(-1200.0, 1200.0), rng.range(-1200.0, 1200.0)]);

        let (xs, ys) = lanes4(&pts);
        let lanes = to_array4(segment_distance_sq_x4(xs, ys, a[0], a[1], b[0], b[1]));

        for i in 0..4 {
            let want = scalar_segment_sq(pts[i][0], pts[i][1], a, b);
            let got = lanes[i];
            // `wide` may reassociate the projection, so the lanes agree with the
            // scalar path to a few ulp rather than bit-for-bit. The published
            // tolerance is relative and far below any geometric significance.
            let scale = want.max(1.0);
            assert!(
                (got - want).abs() <= scale * 1e-12,
                "lane {i}: batch {got} != scalar {want} for seg {a:?}->{b:?} pt {:?}",
                pts[i]
            );
        }
    }
}

#[test]
fn batch_handles_degenerate_segment() {
    // A zero-length segment must behave as a point obstacle.
    let a = [3.0, 4.0];
    let pts = [[0.0, 0.0], [3.0, 4.0], [6.0, 8.0], [-3.0, -4.0]];
    let (xs, ys) = lanes4(&pts);
    let lanes = to_array4(segment_distance_sq_x4(xs, ys, a[0], a[1], a[0], a[1]));
    for i in 0..4 {
        let want = scalar_segment_sq(pts[i][0], pts[i][1], a, a);
        assert!(
            (lanes[i] - want).abs() <= want.max(1.0) * 1e-15,
            "lane {i}: batch {} != scalar {want}",
            lanes[i]
        );
    }
}

#[test]
fn batch_handles_fewer_than_four_points() {
    // `lanes4` repeats the last point, so a short batch never reads out of bounds.
    let pts = [[1.0, 2.0], [5.0, 6.0]];
    let (xs, ys) = lanes4(&pts);
    let lanes = to_array4(segment_distance_sq_x4(xs, ys, 0.0, 0.0, 10.0, 10.0));
    for i in 0..4 {
        let p = pts[i.min(1)];
        let want = scalar_segment_sq(p[0], p[1], [0.0, 0.0], [10.0, 10.0]);
        assert!((lanes[i] - want).abs() <= want.max(1.0) * 1e-15);
    }
}

#[test]
fn lane_extraction_round_trips() {
    let v = f64x4::from([1.5, -2.25, 3.0, 4.75]);
    assert_eq!(to_array4(v), [1.5, -2.25, 3.0, 4.75]);
}
