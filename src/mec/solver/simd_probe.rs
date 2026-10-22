//! Four-lane clearance batch for the climb, behind the `simd` feature.
//!
//! The compass climb probes up to 32 candidate points per round and takes the
//! best. Those probes are independent, so four of them evaluate as one vector
//! operation: the lane arithmetic is identical, only the coordinate differs.
//!
//! The batch computes the distance from four points to one segment. It is used
//! beneath the same BVH descent as the scalar path, so a batched call returns the
//! same numbers the scalar loop would, and the feature is a widening of the inner
//! arithmetic rather than a different algorithm.
//!
//! # Why a segment batch
//!
//! The clearance loop is a reduction over primitives, so vectorising across
//! primitives needs a SoA layout the index does not have. Vectorising across
//! *probes* needs nothing: the four points already exist, and each lane walks the
//! same primitives. That is why the batch is over points, not obstacles.
#![cfg(feature = "simd")]

use wide::f64x4;

/// Squared distance from four points to one segment, lane-wise.
///
/// Returns the squared distance for each point in lane order, matching the scalar
/// `point_segment_distance_sq` exactly: same clamp, same projection, same
/// fallback for a degenerate segment.
#[inline]
pub fn segment_distance_sq_x4(xs: f64x4, ys: f64x4, ax: f64, ay: f64, bx: f64, by: f64) -> f64x4 {
    let ex = bx - ax;
    let ey = by - ay;
    let len_sq = ex * ex + ey * ey;

    if len_sq <= 0.0 {
        // A zero-length segment is a point: the distance is to that point.
        let dx = xs - f64x4::splat(ax);
        let dy = ys - f64x4::splat(ay);
        return dx * dx + dy * dy;
    }

    let dx = xs - f64x4::splat(ax);
    let dy = ys - f64x4::splat(ay);
    let t = (dx * f64x4::splat(ex) + dy * f64x4::splat(ey)) / f64x4::splat(len_sq);
    let t = t.max(f64x4::splat(0.0)).min(f64x4::splat(1.0));

    let px = dx - t * f64x4::splat(ex);
    let py = dy - t * f64x4::splat(ey);
    px * px + py * py
}

/// Fill four lanes from up to four points, repeating the last when short.
#[inline]
pub fn lanes4(points: &[[f64; 2]]) -> (f64x4, f64x4) {
    let mut xs = [0.0f64; 4];
    let mut ys = [0.0f64; 4];
    for i in 0..4 {
        let p = points[i.min(points.len().saturating_sub(1))];
        xs[i] = p[0];
        ys[i] = p[1];
    }
    (f64x4::from(xs), f64x4::from(ys))
}

/// Extract the lanes of a result as a fixed array.
#[inline]
pub fn to_array4(v: f64x4) -> [f64; 4] {
    v.to_array()
}
