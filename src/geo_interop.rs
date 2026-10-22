//! georust ecosystem interop (feature `geo-traits`).
//!
//! Accepts any foreign geometry that implements `GeometryTrait`, converting
//! through `ToGeoGeometry` so callers holding wkb views, geoarrow arrays or
//! geozero streams need no intermediate materialization. Only the first two
//! dimensions are kept.
//!
//! A foreign geometry with no representable conversion, such as an empty point,
//! is rejected with an error rather than panicking.

use geo_traits::to_geo::{ToGeoGeometry, ToGeoPolygon};
use geo_traits::{GeometryTrait, PolygonTrait};
use geo_types::Geometry;

use super::{largest_empty_circle, MecError, MecOptions, MecResult};

/// Solve the largest empty circle on any foreign geometry.
pub fn largest_empty_circle_trait<G>(
    obstacles: &G,
    opts: &MecOptions,
) -> Result<MecResult, MecError>
where
    G: GeometryTrait<T = f64>,
{
    let geom = to_geometry(obstacles)?;
    largest_empty_circle(&geom, None, opts)
}

/// Solve the largest empty circle with a foreign polygonal boundary.
pub fn largest_empty_circle_trait_with_boundary<G, B>(
    obstacles: &G,
    boundary: &B,
    opts: &MecOptions,
) -> Result<MecResult, MecError>
where
    G: GeometryTrait<T = f64>,
    B: PolygonTrait<T = f64>,
{
    let geom = to_geometry(obstacles)?;
    let poly: geo_types::Polygon<f64> = boundary.to_polygon();
    largest_empty_circle(&geom, Some(&Geometry::Polygon(poly)), opts)
}

fn to_geometry<G>(g: &G) -> Result<Geometry<f64>, MecError>
where
    G: GeometryTrait<T = f64>,
{
    g.try_to_geometry().ok_or_else(|| {
        MecError::InvalidInput(
            "foreign geometry has no geo-types equivalent, such as an empty point".to_string(),
        )
    })
}
