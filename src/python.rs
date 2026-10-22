//! PyO3 bindings (feature `python`).
//!
//! Two entry points: point and segment obstacles, with an optional boundary
//! ring. Rings arrive as `(x, y)` tuples. A rejected input, and a solver panic,
//! both raise `ValueError` rather than aborting the interpreter.

use geo_types::{Coord, Geometry, LineString, MultiLineString, MultiPoint, Point, Polygon};
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;

use super::{largest_empty_circle, MecApproach, MecOptions, MecResult};

/// Result of a largest empty circle solve.
#[pyclass(name = "MecResult")]
#[derive(Clone)]
pub struct MecResultPy {
    /// Circle centre x coordinate.
    #[pyo3(get)]
    pub center_x: f64,
    /// Circle centre y coordinate.
    #[pyo3(get)]
    pub center_y: f64,
    /// Circle radius.
    #[pyo3(get)]
    pub radius: f64,
    /// Squared radius.
    #[pyo3(get)]
    pub radius_sq: f64,
    /// "exact_voronoi" or "grid".
    #[pyo3(get)]
    pub approach: String,
    /// Candidates the solve evaluated.
    #[pyo3(get)]
    pub candidate_count: usize,
    /// Obstacle ordinals the circle touches.
    #[pyo3(get)]
    pub support: Vec<usize>,
}

impl From<MecResult> for MecResultPy {
    fn from(r: MecResult) -> Self {
        Self {
            center_x: r.center.x(),
            center_y: r.center.y(),
            radius: r.radius,
            radius_sq: r.radius_sq,
            approach: match r.approach {
                MecApproach::ExactVoronoi => "exact_voronoi".to_string(),
                MecApproach::GridSolver => "grid".to_string(),
            },
            candidate_count: r.candidate_count,
            support: r.support,
        }
    }
}

fn ring_from(pairs: Vec<(f64, f64)>) -> LineString<f64> {
    LineString::from(
        pairs
            .into_iter()
            .map(|(x, y)| Coord { x, y })
            .collect::<Vec<_>>(),
    )
}

fn obstacles_from(
    points: Vec<(f64, f64)>,
    segments: Vec<(f64, f64, f64, f64)>,
) -> PyResult<Geometry<f64>> {
    let mp = MultiPoint::new(points.into_iter().map(|(x, y)| Point::new(x, y)).collect());
    let mls = MultiLineString::new(
        segments
            .into_iter()
            .map(|(x0, y0, x1, y1)| {
                LineString::from(vec![Coord { x: x0, y: y0 }, Coord { x: x1, y: y1 }])
            })
            .collect(),
    );

    let has_points = !mp.0.is_empty();
    let has_segments = !mls.0.is_empty();

    match (has_points, has_segments) {
        (true, false) => Ok(Geometry::MultiPoint(mp)),
        (false, true) => Ok(Geometry::MultiLineString(mls)),
        (true, true) => Ok(Geometry::GeometryCollection(
            geo_types::GeometryCollection::new_from(vec![
                Geometry::MultiPoint(mp),
                Geometry::MultiLineString(mls),
            ]),
        )),
        (false, false) => Err(PyValueError::new_err("no obstacles supplied")),
    }
}

fn solve(
    points: Vec<(f64, f64)>,
    segments: Vec<(f64, f64, f64, f64)>,
    boundary: Option<Vec<(f64, f64)>>,
) -> PyResult<MecResultPy> {
    let geom = obstacles_from(points, segments)?;
    let bgeom = boundary.map(|ring| Geometry::Polygon(Polygon::new(ring_from(ring), Vec::new())));
    match largest_empty_circle(&geom, bgeom.as_ref(), &MecOptions::default()) {
        Ok(r) => Ok(MecResultPy::from(r)),
        Err(e) => Err(PyValueError::new_err(e.to_string())),
    }
}

/// Solve the largest empty circle among points and segments.
///
/// `points` is `(x, y)` tuples. `segments` is `(x0, y0, x1, y1)` tuples.
/// `boundary` is an optional closed ring of `(x, y)` tuples.
#[pyfunction]
#[pyo3(signature = (points, segments=None, boundary=None))]
pub fn solve_mec_py(
    py: Python<'_>,
    points: Vec<(f64, f64)>,
    segments: Option<Vec<(f64, f64, f64, f64)>>,
    boundary: Option<Vec<(f64, f64)>>,
) -> PyResult<MecResultPy> {
    py.allow_threads(|| solve(points, segments.unwrap_or_default(), boundary))
}

/// Register the module contents.
pub fn register(_py: Python<'_>, m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<MecResultPy>()?;
    m.add_function(wrap_pyfunction!(solve_mec_py, m)?)?;
    Ok(())
}
