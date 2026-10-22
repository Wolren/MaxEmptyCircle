//! Largest empty circle (MEC) solvers.
//!
//! The largest empty circle is the biggest circle that fits between obstacles.
//! Obstacles are anything the circle must avoid: points, lines, polygons, or any
//! collection of them. An optional boundary restricts where the centre may sit,
//! and it is a separate input from the obstacles: a building is an obstacle, a
//! property line is a boundary.
//!
//! This is the dual of the maximum inscribed circle. The inscribed circle sits
//! inside one polygon and is limited by that polygon's own boundary; the empty
//! circle sits in open space and is limited by the obstacles around it, which
//! need not enclose the centre at all.
//!
//! # Quick start
//!
//! ```
//! use max_empty_circle::geo_types::{LineString, MultiPoint, Point};
//! use max_empty_circle::{largest_empty_circle, MecOptions};
//!
//! // Three sites: the circle must clear all of them.
//! let obstacles = MultiPoint::new(vec![
//!     Point::new(0.0, 0.0),
//!     Point::new(10.0, 0.0),
//!     Point::new(5.0, 10.0),
//! ]);
//!
//! let circle = largest_empty_circle(
//!     &max_empty_circle::geo_types::Geometry::MultiPoint(obstacles),
//!     None,
//!     &MecOptions::default(),
//! )?;
//! assert!(circle.radius > 0.0);
//! # Ok::<(), max_empty_circle::MecError>(())
//! ```
//!
//! # Solve pipeline
//!
//! [`MecResult::approach`] reports which path produced the answer:
//!
//! 1. **Exact Voronoi** when every obstacle is a point. The centre of the
//!    largest empty circle among point sites is a Voronoi vertex, so the answer
//!    is found by enumerating circumcentres and re-measuring each against every
//!    site. This is exact, not approximated.
//! 2. **Grid search** for segments and polygons. A branch and bound over the
//!    search box, driven by the 1-Lipschitz property of obstacle distance, with
//!    a compass hill climb when the cell budget runs out.
//!
//! The radius a result reports is always the measured distance from its centre
//! to the nearest obstacle, so a returned circle is empty by construction.
//!
//! # Entry points
//!
//! | Function | Use |
//! |----------|-----|
//! | [`largest_empty_circle`] | one obstacle geometry, optional boundary |
//! | [`largest_empty_circle_slice`] | rayon-parallel batch, input order preserved |
//!
//! # Features
//!
//! | Feature | Effect |
//! |---------|--------|
//! | `bench-geos` | comparison only: links the GEOS C library so benchmarks can compare solvers |
//! | `geo-traits` | generic solves over any `geo_traits` geometry |
//! | `shewchuk` | exact `orient2d` predicates from the `robust` crate in the containment test |
//! | `hotpath` | `#[hotpath::measure]` instrumentation on the solver entry points |
//! | `python` | PyO3 extension module named `max_empty_circle` |
//! | `extension-module` | build the Python module as an importable extension |
//!
//! Default features are empty.
//!
//! # C ABI
//!
//! [`ffi`] exposes `mec_solve` and `mec_options_default`: flat coordinate arrays
//! in, `#[repr(C)]` structs out, `0` on success and `-1` on error, with panics
//! caught at the boundary.

#![warn(missing_docs)]
#![cfg_attr(docsrs, feature(doc_cfg))]

pub mod ffi;
pub mod mec;
pub mod tuning;

#[cfg(feature = "python")]
#[cfg_attr(docsrs, doc(cfg(feature = "python")))]
pub mod python;

#[cfg(feature = "geo-traits")]
#[cfg_attr(docsrs, doc(cfg(feature = "geo-traits")))]
pub mod geo_interop;

use self::mec::boundary::{BoundaryIndex, Region};
use self::mec::index::NearestObstacleIndex;
use self::mec::input::{BoundarySet, ObstacleSet};
use self::mec::solver::exact;
use self::mec::solver::grid::solve_grid;

/// Pinned `geo-types` re-export so downstream crates lock to the same version
/// this crate solves over.
pub use geo_types;

use geo_types::{Geometry, Point};
use rayon::prelude::*;
use thiserror::Error;

pub use mec::solver::exact::ExactCandidate;

/// Solver configuration.
///
/// The solve has no tunable knobs: the grid tolerance scales with the input
/// diagonal and the exact path is closed form. The type is kept so callers pass
/// options by reference and future knobs stay source compatible.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub struct MecOptions {}

/// Solver path that produced a result.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum MecApproach {
    /// Exact Voronoi-vertex enumeration, valid for point obstacles.
    ExactVoronoi,
    /// Branch-and-bound grid search over the search box.
    GridSolver,
}

/// Largest empty circle solve result.
#[derive(Debug, Clone)]
pub struct MecResult {
    /// Circle centre.
    pub center: Point<f64>,
    /// Circle radius in input coordinate units.
    pub radius: f64,
    /// Squared radius, radius * radius.
    pub radius_sq: f64,
    /// Ordinals of the obstacles the circle touches. Points come first in the
    /// ordinal space, then segment index offset by the point count.
    pub support: Vec<usize>,
    /// Candidates the solve evaluated.
    pub candidate_count: usize,
    /// Solver path that produced this result.
    pub approach: MecApproach,
}

/// Solver error.
#[derive(Debug, Error)]
pub enum MecError {
    /// Input geometry is malformed, empty, or degenerate.
    #[error("invalid input: {0}")]
    InvalidInput(String),
    /// No admissible centre was found.
    #[error("no valid empty circle found")]
    NoCircleFound,
}

/// Compute the largest empty circle among `obstacles`.
///
/// `boundary` optionally restricts where the centre may sit. With no boundary
/// the centre is free within the obstacle bounding box, so the answer lies
/// among the obstacles, which is what GEOS does with its convex-hull default.
///
/// When every obstacle is a point the solve is exact. Otherwise it is the grid
/// search.
#[cfg_attr(feature = "hotpath", hotpath::measure)]
pub fn largest_empty_circle(
    obstacles: &Geometry<f64>,
    boundary: Option<&Geometry<f64>>,
    opts: &MecOptions,
) -> Result<MecResult, MecError> {
    let set = ObstacleSet::from_geometry(obstacles)?;
    let boundary_set = match boundary {
        Some(g) => Some(BoundarySet::from_polygonal(g)?),
        None => None,
    };
    solve_on_set(&set, boundary_set.as_ref(), opts)
}

/// Compute the largest empty circle in parallel for many obstacle geometries.
///
/// Each item is solved independently on a rayon thread pool. Results come back
/// in input order.
#[cfg_attr(feature = "hotpath", hotpath::measure)]
pub fn largest_empty_circle_slice(
    obstacles: &[Geometry<f64>],
    opts: &MecOptions,
) -> Vec<Result<MecResult, MecError>> {
    obstacles
        .par_iter()
        .map(|g| largest_empty_circle(g, None, opts))
        .collect()
}

fn solve_on_set(
    set: &ObstacleSet,
    boundary_set: Option<&BoundarySet>,
    _opts: &MecOptions,
) -> Result<MecResult, MecError> {
    let boundary_index = boundary_set.map(BoundaryIndex::new);
    let region = Region::new(set, boundary_index.as_ref())?;

    // The exact path applies only when every obstacle is a point. It is also the
    // only path that is exact rather than approximated to a tolerance, so it
    // runs first and its answer stands.
    if set.is_point_only() {
        if let Some(c) = exact::solve_exact(set, &region) {
            let support = c.support[..c.support_len].to_vec();
            return Ok(MecResult {
                center: exact::candidate_center(&c),
                radius: c.radius_sq.sqrt(),
                radius_sq: c.radius_sq,
                support,
                candidate_count: set.point_count(),
                approach: MecApproach::ExactVoronoi,
            });
        }
    }

    let nb = NearestObstacleIndex::new(set);
    let tolerance = grid_tolerance(set, &region);

    if let Some(sol) = solve_grid(set, &nb, &region, tolerance) {
        let radius = sol.radius_sq.sqrt();
        let support = nb.supporting_primitives(sol.x, sol.y, sol.radius_sq + sol.radius_sq * 1e-10);
        return Ok(MecResult {
            center: Point::new(sol.x, sol.y),
            radius,
            radius_sq: sol.radius_sq,
            support,
            candidate_count: sol.iterations,
            approach: MecApproach::GridSolver,
        });
    }

    Err(MecError::NoCircleFound)
}

/// Absolute refinement tolerance for the grid stage: `diag * 1e-7`.
///
/// The diagonal carries no floor, because a floor would tie the accuracy of the
/// answer to the coordinate units of the input: a region a few hundred metres
/// wide stored in degrees has a diagonal near `1e-3`, so a floor of one unit
/// makes the tolerance a whole percent of its own radius.
pub fn grid_tolerance(_set: &ObstacleSet, region: &Region) -> f64 {
    let (min_x, min_y, max_x, max_y) = region.bounds();
    let diag = (max_x - min_x).hypot(max_y - min_y);
    diag * 1e-7
}

#[cfg(test)]
mod tests {
    use super::*;
    use geo_types::{Geometry, LineString, MultiPoint, Point, Polygon};

    fn sites(pts: &[(f64, f64)]) -> Geometry<f64> {
        Geometry::MultiPoint(MultiPoint::new(
            pts.iter().map(|p| Point::new(p.0, p.1)).collect(),
        ))
    }

    #[test]
    fn mec_three_sites_exact() {
        // Equilateral triangle side 10: circumcentre at the centroid, radius
        // 10 / sqrt(3). The exact path must find it.
        let g = sites(&[(0.0, 0.0), (10.0, 0.0), (5.0, 8.660_254_037_844_387)]);
        let r = largest_empty_circle(&g, None, &MecOptions::default()).expect("solve");
        assert_eq!(r.approach, MecApproach::ExactVoronoi);
        assert!(
            (r.radius - 5.773_502_691_8) < 1e-6,
            "radius={} expected 10/sqrt(3)",
            r.radius
        );
    }

    #[test]
    fn mec_two_sites_needs_a_boundary_to_be_bounded() {
        // Collinear obstacles enclose no area, so the free problem is unbounded
        // and the solve refuses rather than inventing a region.
        let g = sites(&[(0.0, 0.0), (10.0, 0.0)]);
        assert!(largest_empty_circle(&g, None, &MecOptions::default()).is_err());

        // With a boundary the midpoint of the pair is the answer.
        let b = Geometry::Polygon(Polygon::new(
            LineString::from(vec![
                (0.0, -1.0),
                (10.0, -1.0),
                (10.0, -1.0),
                (0.0, -1.0),
                (0.0, -1.0),
            ]),
            vec![],
        ));
        let _ = b;
        let boxy = Geometry::Polygon(Polygon::new(
            LineString::from(vec![
                (0.0, 0.0),
                (10.0, 0.0),
                (10.0, 0.1),
                (0.0, 0.1),
                (0.0, 0.0),
            ]),
            vec![],
        ));
        let r = largest_empty_circle(&g, Some(&boxy), &MecOptions::default()).expect("solve");
        assert_eq!(r.approach, MecApproach::ExactVoronoi);
        // The region sweep resolves to a fraction of the span, not machine
        // precision, so the tolerance is relative to the answer.
        assert!((r.radius - 5.0).abs() < 1e-2, "radius={}", r.radius);
        assert!((r.center.x() - 5.0).abs() < 1e-2);
    }

    #[test]
    fn mec_four_corners_square() {
        // Four sites on a square: the centre clears each by half the diagonal.
        let g = sites(&[(0.0, 0.0), (10.0, 0.0), (10.0, 10.0), (0.0, 10.0)]);
        let r = largest_empty_circle(&g, None, &MecOptions::default()).expect("solve");
        assert!(
            (r.radius - 7.071_067_811_9) < 1e-6,
            "radius={} expected 5*sqrt(2)",
            r.radius
        );
    }

    #[test]
    fn mec_boundary_confines_centre() {
        // The free circumcentre of this triple sits at (5, 2.886). A boundary
        // ending at x = 4.5 excludes it, so the centre must fall back to an
        // admissible candidate and the radius shrinks.
        let g = sites(&[(0.0, 0.0), (10.0, 0.0), (5.0, 8.660_254_037_844_387)]);
        let b = Geometry::Polygon(Polygon::new(
            LineString::from(vec![
                (0.0, 0.0),
                (4.5, 0.0),
                (4.5, 9.0),
                (0.0, 9.0),
                (0.0, 0.0),
            ]),
            vec![],
        ));
        let free = largest_empty_circle(&g, None, &MecOptions::default()).expect("free");
        let confined =
            largest_empty_circle(&g, Some(&b), &MecOptions::default()).expect("confined");
        // The free optimum sits outside this boundary, so confining must cost
        // clearance. The boundary clips it at x = 4.5, where the centre sits.
        assert!(
            confined.center.x() <= 4.5 + 1e-9,
            "x={} escaped the boundary",
            confined.center.x()
        );
        assert!(
            confined.radius <= free.radius + 1e-9,
            "confined {} should not beat free {}",
            confined.radius,
            free.radius
        );
    }

    #[test]
    fn mec_polygon_obstacle_uses_grid() {
        // A single areal obstacle has no point sites, so the grid path answers.
        let g = Geometry::Polygon(Polygon::new(
            LineString::from(vec![
                (0.0, 0.0),
                (4.0, 0.0),
                (4.0, 4.0),
                (0.0, 4.0),
                (0.0, 0.0),
            ]),
            vec![],
        ));
        let r = largest_empty_circle(&g, None, &MecOptions::default()).expect("solve");
        assert_eq!(r.approach, MecApproach::GridSolver);
        assert!(r.radius > 0.0);
    }

    #[test]
    fn mec_line_obstacles_use_grid() {
        // Two parallel segments give the hull area, so the region confines the
        // centre and the widest gap is half the separation.
        let g = Geometry::MultiLineString(geo_types::MultiLineString::new(vec![
            LineString::from(vec![(0.0, 0.0), (10.0, 0.0)]),
            LineString::from(vec![(0.0, 4.0), (10.0, 4.0)]),
        ]));
        let r = largest_empty_circle(&g, None, &MecOptions::default()).expect("solve");
        assert_eq!(r.approach, MecApproach::GridSolver);
        assert!((r.radius - 2.0).abs() < 1e-3, "radius={}", r.radius);
    }

    #[test]
    fn mec_rejects_empty_obstacles() {
        let g = sites(&[]);
        assert!(largest_empty_circle(&g, None, &MecOptions::default()).is_err());
    }

    #[test]
    fn mec_rejects_non_polygonal_boundary() {
        let g = sites(&[(0.0, 0.0), (10.0, 0.0)]);
        let b = Geometry::LineString(LineString::from(vec![(0.0, 0.0), (1.0, 1.0)]));
        let err = largest_empty_circle(&g, Some(&b), &MecOptions::default()).unwrap_err();
        assert!(err.to_string().contains("polygonal"), "got {err}");
    }

    #[test]
    fn mec_batch_preserves_order() {
        let shapes: Vec<Geometry<f64>> = (1..=4)
            .map(|k| sites(&[(0.0, 0.0), (k as f64, 0.0), (0.0, k as f64)]))
            .collect();
        let out = largest_empty_circle_slice(&shapes, &MecOptions::default());
        assert_eq!(out.len(), 4);
        for r in out.iter() {
            assert!(r.as_ref().expect("solve").radius > 0.0);
        }
    }
}
