//! GEOS comparison through hand-declared C bindings, behind `bench-geos`.
//!
//! The safe `geos` crate never wrapped `GEOSLargestEmptyCircle`, and both it and
//! `geos-sys` pull conversion crates or conflicting binding features, so this
//! module declares the seven C functions it needs and links `geos_c` directly.
//! That keeps the comparison to one dependency-free module and lets it call the
//! exact function being measured.
//!
//! It is never part of the solve pipeline: nothing in [`crate::mec`] calls into
//! it, and the default feature set does not compile it.
//!
//! GEOS runs the same cell-refinement branch and bound as
//! [`crate::mec::solver::grid`], so agreement means two implementations of one
//! algorithm match, not that an independent method confirmed the answer. For
//! point obstacles the native exact path *is* independent of GEOS, so there the
//! comparison does test something: a Voronoi enumeration against an iterative
//! approximation.
//!
//! GEOS defaults the boundary to the convex hull of the obstacles, which is the
//! same default this crate applies, so the two search the same region.
//!
//! # Linking
//!
//! The build script needs `GEOS_LIB_DIR` and `GEOS_VERSION` unless `geos-config`
//! is on `PATH`. On a Miniconda install that is
//! `GEOS_LIB_DIR=D:\Miniconda\Library\lib` and `GEOS_VERSION=3.14.1`.
//!
//! # Safety
//!
//! Every raw pointer is wrapped in one place and freed by a guard. Failure at any
//! step returns an error rather than a partial result.

use geo_types::{Coord, Geometry, LineString, MultiLineString, MultiPoint, Point, Polygon};
use std::ffi::{c_char, c_double, c_int, c_uint, c_void, CString};
use std::ptr;

use super::super::input::ObstacleSet;
use crate::{MecError, MecResult};

/// Opaque GEOS geometry.
#[repr(C)]
struct GEOSGeometry {
    _private: [u8; 0],
}

/// Opaque GEOS coordinate sequence.
#[repr(C)]
struct GEOSCoordSequence {
    _private: [u8; 0],
}

/// Opaque GEOS context handle.
#[repr(C)]
struct GEOSContextHandle {
    _private: [u8; 0],
}

type Context = *mut GEOSContextHandle;

extern "C" {
    fn GEOS_init_r() -> Context;
    fn GEOS_finish_r(handle: Context);
    fn GEOSWKTReader_create_r(handle: Context) -> *mut c_void;
    fn GEOSWKTReader_read_r(
        handle: Context,
        reader: *mut c_void,
        wkt: *const c_char,
    ) -> *mut GEOSGeometry;
    fn GEOSWKTReader_destroy_r(handle: Context, reader: *mut c_void);
    fn GEOSGeom_destroy_r(handle: Context, geom: *mut GEOSGeometry);
    fn GEOSGeom_getCoordSeq_r(
        handle: Context,
        geom: *const GEOSGeometry,
    ) -> *const GEOSCoordSequence;
    fn GEOSCoordSeq_getSize_r(
        handle: Context,
        seq: *const GEOSCoordSequence,
        size: *mut c_uint,
    ) -> c_int;
    fn GEOSCoordSeq_getX_r(
        handle: Context,
        seq: *const GEOSCoordSequence,
        idx: c_uint,
        val: *mut c_double,
    ) -> c_int;
    fn GEOSCoordSeq_getY_r(
        handle: Context,
        seq: *const GEOSCoordSequence,
        idx: c_uint,
        val: *mut c_double,
    ) -> c_int;
    fn GEOSLargestEmptyCircle_r(
        handle: Context,
        obstacles: *const GEOSGeometry,
        boundary: *const GEOSGeometry,
        tolerance: c_double,
    ) -> *mut GEOSGeometry;
}

/// Tolerance to hand GEOS for an obstacle set, scaled by its diagonal.
///
/// The native grid refines to `diag * 1e-7`, so this hands GEOS `diag * 1e-8`,
/// floored at `1e-12`: a tighter GEOS is given the work, and the comparison then
/// measures the solvers rather than the tolerance gap between them.
pub fn geos_tolerance_for(set: &ObstacleSet) -> f64 {
    let Some((min_x, min_y, max_x, max_y)) = set.bounds() else {
        return 1e-6;
    };
    let diag = (max_x - min_x).hypot(max_y - min_y).max(1.0);
    (diag * 1e-8).max(1e-12)
}

/// The centre and radius GEOS reports, with the radius re-measured natively.
#[derive(Debug, Clone, Copy)]
pub struct GeosCircle {
    /// Centre x.
    pub x: f64,
    /// Centre y.
    pub y: f64,
    /// Radius, measured against the obstacle set.
    pub radius: f64,
}

/// A GEOS geometry that frees itself.
struct OwnedGeom {
    ptr: *mut GEOSGeometry,
    handle: Context,
}

impl Drop for OwnedGeom {
    fn drop(&mut self) {
        if !self.ptr.is_null() && !self.handle.is_null() {
            unsafe { GEOSGeom_destroy_r(self.handle, self.ptr) };
        }
    }
}

/// A GEOS context that finishes itself, created per call.
///
/// Per call rather than shared, because the handle is not thread safe and the
/// comparison harness solves on several threads.
struct OwnedContext(Context);

impl OwnedContext {
    fn new() -> Result<Self, MecError> {
        let h = unsafe { GEOS_init_r() };
        if h.is_null() {
            return Err(MecError::GeosFailed("GEOS_init_r failed".to_string()));
        }
        Ok(Self(h))
    }

    fn get(&self) -> Context {
        self.0
    }
}

impl Drop for OwnedContext {
    fn drop(&mut self) {
        if !self.0.is_null() {
            unsafe { GEOS_finish_r(self.0) };
        }
    }
}

/// Build a WKT string for any geometry.
///
/// WKT rather than the C geometry constructors: one code path for every variant,
/// and the string is what a caller would hand GEOS anyway.
fn to_wkt(geom: &Geometry<f64>) -> String {
    match geom {
        Geometry::Point(p) => point_wkt(p, "POINT "),
        Geometry::MultiPoint(mp) => {
            let parts: Vec<String> = mp.0.iter().map(|p| point_wkt(p, "")).collect();
            format!("MULTIPOINT ({})", parts.join(", "))
        }
        Geometry::Line(l) => format!(
            "LINESTRING ({} {}, {} {})",
            l.start.x, l.start.y, l.end.x, l.end.y
        ),
        Geometry::LineString(ls) => line_string_wkt(ls, "LINESTRING "),
        Geometry::MultiLineString(mls) => {
            let parts: Vec<String> = mls
                .0
                .iter()
                .map(|ls| line_string_wkt(ls, "").trim().to_string())
                .collect();
            format!("MULTILINESTRING ({})", parts.join(", "))
        }
        Geometry::Polygon(p) => polygon_wkt(p, "POLYGON "),
        Geometry::MultiPolygon(mp) => {
            let parts: Vec<String> =
                mp.0.iter()
                    .map(|p| polygon_wkt(p, "").trim().to_string())
                    .collect();
            format!("MULTIPOLYGON ({})", parts.join(", "))
        }
        Geometry::Rect(r) => polygon_wkt(&r.to_polygon(), "POLYGON "),
        Geometry::Triangle(t) => polygon_wkt(&t.to_polygon(), "POLYGON "),
        Geometry::GeometryCollection(gc) => {
            let parts: Vec<String> = gc.0.iter().map(to_wkt).collect();
            format!("GEOMETRYCOLLECTION ({})", parts.join(", "))
        }
    }
}

fn point_wkt(p: &Point<f64>, prefix: &str) -> String {
    format!("{prefix}({} {})", p.x(), p.y())
}

fn line_string_wkt(ls: &LineString<f64>, prefix: &str) -> String {
    let coords: Vec<String> = ls.0.iter().map(|c| format!("{} {}", c.x, c.y)).collect();
    format!("{prefix}({})", coords.join(", "))
}

fn polygon_wkt(p: &Polygon<f64>, prefix: &str) -> String {
    let mut rings: Vec<String> = Vec::with_capacity(1 + p.interiors().len());
    rings.push(ring_body(p.exterior()));
    for hole in p.interiors() {
        rings.push(ring_body(hole));
    }
    format!("{prefix}({})", rings.join(", "))
}

fn ring_body(ls: &LineString<f64>) -> String {
    let coords: Vec<String> = ls.0.iter().map(|c| format!("{} {}", c.x, c.y)).collect();
    format!("({})", coords.join(", "))
}

/// Parse a WKT geometry into a self-freeing handle.
fn parse_wkt(h: Context, wkt: &str) -> Result<OwnedGeom, MecError> {
    let c = CString::new(wkt)
        .map_err(|_| MecError::GeosFailed("WKT contains a NUL byte".to_string()))?;
    let reader = unsafe { GEOSWKTReader_create_r(h) };
    if reader.is_null() {
        return Err(MecError::GeosFailed(
            "GEOSWKTReader_create_r failed".to_string(),
        ));
    }
    let geom = unsafe { GEOSWKTReader_read_r(h, reader, c.as_ptr()) };
    unsafe { GEOSWKTReader_destroy_r(h, reader) };
    if geom.is_null() {
        return Err(MecError::GeosFailed(format!(
            "GEOS failed to parse WKT: {wkt}"
        )));
    }
    Ok(OwnedGeom {
        ptr: geom,
        handle: h,
    })
}

/// Read the first coordinate of a geometry as (x, y).
fn first_coord(h: Context, geom: *const GEOSGeometry) -> Option<(f64, f64)> {
    unsafe {
        let seq = GEOSGeom_getCoordSeq_r(h, geom);
        if seq.is_null() {
            return None;
        }
        let mut size: c_uint = 0;
        if GEOSCoordSeq_getSize_r(h, seq, &mut size) == 0 || size == 0 {
            return None;
        }
        let mut x = 0.0;
        let mut y = 0.0;
        if GEOSCoordSeq_getX_r(h, seq, 0, &mut x) == 0 {
            return None;
        }
        if GEOSCoordSeq_getY_r(h, seq, 0, &mut y) == 0 {
            return None;
        }
        Some((x, y))
    }
}

/// Solve with GEOS at an explicit tolerance.
pub fn solve_geos_with_tolerance(
    set: &ObstacleSet,
    boundary: Option<&Geometry<f64>>,
    tolerance: f64,
) -> Result<GeosCircle, MecError> {
    let ctx = OwnedContext::new()?;
    let h = ctx.get();

    let obstacles_wkt = obstacles_to_wkt(set)?;
    let obstacles = parse_wkt(h, &obstacles_wkt)?;

    let boundary_owned = match boundary {
        Some(b) => Some(parse_wkt(h, &to_wkt(b))?),
        None => None,
    };
    let boundary_ptr = boundary_owned
        .as_ref()
        .map(|g| g.ptr as *const GEOSGeometry)
        .unwrap_or(ptr::null());

    let line = unsafe { GEOSLargestEmptyCircle_r(h, obstacles.ptr, boundary_ptr, tolerance) };
    if line.is_null() {
        return Err(MecError::GeosFailed(
            "GEOSLargestEmptyCircle_r returned null".to_string(),
        ));
    }
    let line = OwnedGeom {
        ptr: line,
        handle: h,
    };

    let (x, y) = first_coord(h, line.ptr)
        .ok_or_else(|| MecError::GeosFailed("could not read GEOS circle centre".to_string()))?;

    // Measure against the native table rather than trusting the radius line, so
    // both solvers report radii on the same basis.
    let radius = set.clearance_sq(x, y).max(0.0).sqrt();
    Ok(GeosCircle { x, y, radius })
}

/// Solve with GEOS at the default comparison tolerance.
pub fn solve_geos(
    set: &ObstacleSet,
    boundary: Option<&Geometry<f64>>,
) -> Result<GeosCircle, MecError> {
    let tolerance = geos_tolerance_for(set);
    solve_geos_with_tolerance(set, boundary, tolerance)
}

/// Build the WKT for an obstacle set: points as multipoint, edges as multilinestring.
fn obstacles_to_wkt(set: &ObstacleSet) -> Result<String, MecError> {
    let mp = MultiPoint::new(
        set.points
            .iter()
            .map(|p| Point::new(p[0], p[1]))
            .collect::<Vec<_>>(),
    );
    let mut segs: Vec<LineString<f64>> = Vec::with_capacity(set.segments.len());
    for i in 0..set.segments.len() {
        let s = &set.segments;
        segs.push(LineString::from(vec![
            Coord {
                x: s.ax[i],
                y: s.ay[i],
            },
            Coord {
                x: s.bx[i],
                y: s.by[i],
            },
        ]));
    }
    let mls = MultiLineString::new(segs);

    let has_points = !mp.0.is_empty();
    let has_segments = !mls.0.is_empty();

    match (has_points, has_segments) {
        (true, false) => Ok(to_wkt(&Geometry::MultiPoint(mp))),
        (false, true) => Ok(to_wkt(&Geometry::MultiLineString(mls))),
        (true, true) => Ok(to_wkt(&Geometry::GeometryCollection(
            geo_types::GeometryCollection::new_from(vec![
                Geometry::MultiPoint(mp),
                Geometry::MultiLineString(mls),
            ]),
        ))),
        (false, false) => Err(MecError::InvalidInput("obstacle set is empty".to_string())),
    }
}

/// Wrap a GEOS circle as a native result, for comparison harnesses.
pub fn geos_circle_to_result(c: GeosCircle) -> MecResult {
    MecResult {
        center: Point::new(c.x, c.y),
        radius: c.radius,
        radius_sq: c.radius * c.radius,
        support: Vec::new(),
        candidate_count: 1,
        approach: crate::MecApproach::Geos,
    }
}

/// A GEOS context with the obstacle geometry already parsed.
///
/// Build once, then time [`PreparedGeos::solve`] alone. The plain
/// [`solve_geos_with_tolerance`] builds the WKT string and parses it inside its
/// own call, which makes a repeated timing measure parsing rather than solving.
pub struct PreparedGeos {
    context: OwnedContext,
    obstacles: OwnedGeom,
    boundary: Option<OwnedGeom>,
}

// The GEOS context handle is owned by this struct and never shared, so moving the
// prepared pair to another thread moves the only reference to it.
unsafe impl Send for PreparedGeos {}

impl PreparedGeos {
    /// Build the GEOS geometry for an obstacle set and optional boundary.
    pub fn new(set: &ObstacleSet, boundary: Option<&Geometry<f64>>) -> Result<Self, MecError> {
        let context = OwnedContext::new()?;
        let h = context.get();
        let obstacles = parse_wkt(h, &obstacles_to_wkt(set)?)?;
        let boundary = match boundary {
            Some(b) => Some(parse_wkt(h, &to_wkt(b))?),
            None => None,
        };
        Ok(Self {
            context,
            obstacles,
            boundary,
        })
    }

    /// Solve at an explicit tolerance, with no setup inside the call.
    ///
    /// `set` supplies the measured radius, so both solvers report on one basis.
    pub fn solve(&self, set: &ObstacleSet, tolerance: f64) -> Result<GeosCircle, MecError> {
        let h = self.context.get();
        let boundary_ptr = self
            .boundary
            .as_ref()
            .map(|g| g.ptr as *const GEOSGeometry)
            .unwrap_or(ptr::null());

        let line =
            unsafe { GEOSLargestEmptyCircle_r(h, self.obstacles.ptr, boundary_ptr, tolerance) };
        if line.is_null() {
            return Err(MecError::GeosFailed(
                "GEOSLargestEmptyCircle_r returned null".to_string(),
            ));
        }
        let line = OwnedGeom {
            ptr: line,
            handle: h,
        };

        let (x, y) = first_coord(h, line.ptr)
            .ok_or_else(|| MecError::GeosFailed("could not read GEOS circle centre".to_string()))?;
        let radius = set.clearance_sq(x, y).max(0.0).sqrt();
        Ok(GeosCircle { x, y, radius })
    }
}
