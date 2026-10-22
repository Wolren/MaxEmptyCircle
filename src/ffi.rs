//! C ABI for the MEC solvers.
//!
//! Every exported function returns 0 on success, -1 on error, and never lets a
//! Rust panic cross into the host process. Obstacles arrive as a flat point
//! array plus a flat segment array; the boundary is optional and polygonal.

use super::{largest_empty_circle, MecApproach, MecOptions, MecResult};
use geo_types::{Coord, Geometry, LineString, MultiLineString, MultiPoint, Point, Polygon};
use libc::{c_double, c_int, size_t};
use std::slice;

/// MEC result.
#[repr(C)]
#[derive(Clone, Copy, Default)]
pub struct MecResultC {
    /// Circle centre x coordinate.
    pub center_x: c_double,
    /// Circle centre y coordinate.
    pub center_y: c_double,
    /// Circle radius in input coordinate units.
    pub radius: c_double,
    /// Squared radius, radius * radius.
    pub radius_sq: c_double,
    /// 1 = exact Voronoi, 2 = grid.
    pub approach: c_int,
    /// Candidates the solve evaluated.
    pub candidate_count: size_t,
}

/// MEC options.
///
/// Reserved: the solve has no tunables, so this struct carries no fields today.
/// It is kept so the C ABI has a stable options pointer to pass.
#[repr(C)]
#[derive(Clone, Copy, Default)]
pub struct MecOptionsC {
    _reserved: c_int,
}

/// Decode a flat [x0, y0, x1, y1, ...] array into pairs. Validates length.
///
/// # Safety
/// `coords` must point to `coords_len` readable doubles, or be null
/// (null returns None).
pub unsafe fn decode_coords(
    coords: *const c_double,
    coords_len: size_t,
) -> Option<Vec<(f64, f64)>> {
    if coords.is_null() || coords_len == 0 || !coords_len.is_multiple_of(2) {
        return None;
    }
    let s = slice::from_raw_parts(coords, coords_len);
    let mut out = Vec::with_capacity(coords_len / 2);
    for i in (0..coords_len).step_by(2) {
        let (x, y) = (s[i], s[i + 1]);
        if !x.is_finite() || !y.is_finite() {
            return None;
        }
        out.push((x, y));
    }
    Some(out)
}

/// Run a closure behind the panic boundary. Panic becomes -1.
pub fn guard(f: impl FnOnce() -> c_int + std::panic::UnwindSafe) -> c_int {
    std::panic::catch_unwind(f).unwrap_or(-1)
}

fn approach_to(v: MecApproach) -> c_int {
    match v {
        MecApproach::ExactVoronoi => 1,
        MecApproach::GridSolver => 2,
        MecApproach::Geos => 3,
    }
}

fn fill(r: &MecResult, out: *mut MecResultC) -> c_int {
    if out.is_null() {
        return -1;
    }
    unsafe {
        *out = MecResultC {
            center_x: r.center.x(),
            center_y: r.center.y(),
            radius: r.radius,
            radius_sq: r.radius_sq,
            approach: approach_to(r.approach),
            candidate_count: r.candidate_count as size_t,
        };
    }
    0
}

/// Build the obstacle geometry from the flat C arrays.
///
/// Points go through `MultiPoint` and segments through `MultiLineString`, so
/// the public constructor builds both tables in one place and the C path cannot
/// drift from the Rust one.
unsafe fn build_obstacles(
    points: *const c_double,
    points_len: size_t,
    segments: *const c_double,
    segments_len: size_t,
) -> Option<Geometry<f64>> {
    let points = if points_len > 0 {
        Some(decode_coords(points, points_len)?)
    } else {
        None
    };

    let segment_strings: Vec<LineString<f64>> = if segments_len > 0 {
        if !segments_len.is_multiple_of(4) {
            return None;
        }
        let s = slice::from_raw_parts(segments, segments_len);
        let mut out = Vec::with_capacity(segments_len / 4);
        for w in s.chunks(4) {
            if !w.iter().all(|v| v.is_finite()) {
                return None;
            }
            out.push(LineString::from(vec![
                Coord { x: w[0], y: w[1] },
                Coord { x: w[2], y: w[3] },
            ]));
        }
        out
    } else {
        Vec::new()
    };

    match (points, segment_strings.is_empty()) {
        (Some(pts), true) => Some(Geometry::MultiPoint(MultiPoint::new(
            pts.into_iter().map(|p| Point::new(p.0, p.1)).collect(),
        ))),
        (None, false) => Some(Geometry::MultiLineString(MultiLineString::new(
            segment_strings,
        ))),
        (Some(pts), false) => {
            let parts: Vec<Geometry<f64>> = vec![
                Geometry::MultiLineString(MultiLineString::new(segment_strings)),
                Geometry::MultiPoint(MultiPoint::new(
                    pts.into_iter().map(|p| Point::new(p.0, p.1)).collect(),
                )),
            ];
            Some(Geometry::GeometryCollection(
                geo_types::GeometryCollection::new_from(parts),
            ))
        }
        (None, true) => None,
    }
}

/// Solve the largest empty circle.
///
/// `points` is a flat [x, y, ...] array of point obstacles. `segments` is a flat
/// [x0, y0, x1, y1, ...] array of segment obstacles, four doubles per segment.
/// `boundary` is an optional flat closed ring; pass null and 0 for no boundary.
///
/// Returns 0 on success, -1 on error.
///
/// # Safety
/// All pointer arguments must be valid for reads of their declared lengths for
/// the duration of the call. Null and short inputs return -1 instead of
/// crashing.
#[no_mangle]
pub unsafe extern "C" fn mec_solve(
    points: *const c_double,
    points_len: size_t,
    segments: *const c_double,
    segments_len: size_t,
    boundary: *const c_double,
    boundary_len: size_t,
    options: *const MecOptionsC,
    result: *mut MecResultC,
) -> c_int {
    guard(|| {
        if result.is_null() {
            return -1;
        }

        let Some(geom) = build_obstacles(points, points_len, segments, segments_len) else {
            return -1;
        };

        let bgeom = if boundary_len > 0 {
            let Some(ring) = decode_coords(boundary, boundary_len) else {
                return -1;
            };
            if ring.len() < 4 {
                return -1;
            }
            Some(Geometry::Polygon(Polygon::new(
                LineString::from(
                    ring.into_iter()
                        .map(|p| Coord { x: p.0, y: p.1 })
                        .collect::<Vec<_>>(),
                ),
                vec![],
            )))
        } else {
            None
        };

        // Options carry nothing today; reading the pointer keeps the ABI shape.
        let _ = if options.is_null() {
            MecOptionsC::default()
        } else {
            *options
        };

        match largest_empty_circle(&geom, bgeom.as_ref(), &MecOptions::default()) {
            Ok(r) => fill(&r, result),
            Err(_) => -1,
        }
    })
}

#[no_mangle]
/// Return the default C ABI options.
pub extern "C" fn mec_options_default() -> MecOptionsC {
    MecOptionsC::default()
}
