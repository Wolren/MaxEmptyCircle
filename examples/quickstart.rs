//! Solve the largest empty circle for a set of obstacles, then with a boundary.
//!
//! Run with `cargo run --example quickstart`.

use max_empty_circle::geo_types::{
    Geometry, LineString, MultiLineString, MultiPoint, Point, Polygon,
};
use max_empty_circle::{largest_empty_circle, MecOptions};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let opts = MecOptions::default();

    // Three sites: the exact Voronoi path answers.
    let sites = Geometry::MultiPoint(MultiPoint::new(vec![
        Point::new(0.0, 0.0),
        Point::new(10.0, 0.0),
        Point::new(5.0, 8.660_254_037_844_387),
    ]));
    let circle = largest_empty_circle(&sites, None, &opts)?;
    println!(
        "three sites: centre ({:.4}, {:.4}) radius {:.4} via {:?}",
        circle.center.x(),
        circle.center.y(),
        circle.radius,
        circle.approach
    );

    // Two parallel walls: the grid path answers.
    let walls = Geometry::MultiLineString(MultiLineString::new(vec![
        LineString::from(vec![(0.0, 0.0), (100.0, 0.0)]),
        LineString::from(vec![(0.0, 50.0), (100.0, 50.0)]),
    ]));
    let circle = largest_empty_circle(&walls, None, &opts)?;
    println!(
        "two walls: centre ({:.4}, {:.4}) radius {:.4} via {:?}",
        circle.center.x(),
        circle.center.y(),
        circle.radius,
        circle.approach
    );

    // Confine the centre to a boundary smaller than the obstacle span.
    let boundary = Geometry::Polygon(Polygon::new(
        LineString::from(vec![
            (0.0, 0.0),
            (4.0, 0.0),
            (4.0, 9.0),
            (0.0, 9.0),
            (0.0, 0.0),
        ]),
        vec![],
    ));
    let circle = largest_empty_circle(&sites, Some(&boundary), &opts)?;
    println!(
        "confined:  centre ({:.4}, {:.4}) radius {:.4} via {:?}",
        circle.center.x(),
        circle.center.y(),
        circle.radius,
        circle.approach
    );

    Ok(())
}
