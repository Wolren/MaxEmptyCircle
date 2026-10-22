//! Obstacle and boundary normalization for the MEC solvers.
//!
//! [`ObstacleSet`] flattens every obstacle geometry into two primitive tables:
//! points, which the exact Voronoi path consumes, and segments, which the grid
//! path consumes. [`SegmentIndex`] is the struct-of-arrays edge table carried
//! over from the sibling `max-inscribed-circle` crate, because the distance
//! query against a set of edges is the same query in both problems.
use geo_types::{Coord, Geometry, Line, LineString, MultiPolygon, Point, Polygon};

use crate::MecError;

const DEGENERATE_EPS: f64 = 1e-12;

/// Struct-of-arrays segment table over every obstacle edge.
#[derive(Debug, Clone, Default)]
pub struct SegmentIndex {
    /// Segment start x; index i matches index i of every other array here.
    pub ax: Vec<f64>,
    /// Segment start y, aligned with `ax`.
    pub ay: Vec<f64>,
    /// Segment end x, aligned with `ax`.
    pub bx: Vec<f64>,
    /// Segment end y, aligned with `ax`.
    pub by: Vec<f64>,
    /// Index of the input geometry that owns each segment.
    pub owner: Vec<usize>,
    /// Minimal x of each segment's bounding box.
    pub bbox_minx: Vec<f64>,
    /// Maximal x of each segment's bounding box.
    pub bbox_maxx: Vec<f64>,
    /// Minimal y of each segment's bounding box.
    pub bbox_miny: Vec<f64>,
    /// Maximal y of each segment's bounding box.
    pub bbox_maxy: Vec<f64>,
    /// Segment direction x, computed as bx - ax.
    pub dir_x: Vec<f64>,
    /// Segment direction y, computed as by - ay.
    pub dir_y: Vec<f64>,
    /// Reciprocal squared segment length.
    pub inv_len_sq: Vec<f64>,
}

impl SegmentIndex {
    /// Number of segments in the table.
    pub fn len(&self) -> usize {
        self.ax.len()
    }

    /// True when the table holds no segments.
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    fn push(&mut self, a: [f64; 2], b: [f64; 2], owner: usize) {
        let dx = b[0] - a[0];
        let dy = b[1] - a[1];
        let len_sq = dx * dx + dy * dy;
        if len_sq <= DEGENERATE_EPS {
            return;
        }
        self.ax.push(a[0]);
        self.ay.push(a[1]);
        self.bx.push(b[0]);
        self.by.push(b[1]);
        self.owner.push(owner);
        self.bbox_minx.push(a[0].min(b[0]));
        self.bbox_maxx.push(a[0].max(b[0]));
        self.bbox_miny.push(a[1].min(b[1]));
        self.bbox_maxy.push(a[1].max(b[1]));
        self.dir_x.push(dx);
        self.dir_y.push(dy);
        self.inv_len_sq.push(1.0 / len_sq);
    }

    #[inline]
    /// Squared distance from point (x, y) to segment `seg_idx`.
    pub fn point_segment_distance_sq(&self, seg_idx: usize, x: f64, y: f64) -> f64 {
        let ax = self.ax[seg_idx];
        let ay = self.ay[seg_idx];
        let dx = self.dir_x[seg_idx];
        let dy = self.dir_y[seg_idx];
        let t = (((x - ax) * dx + (y - ay) * dy) * self.inv_len_sq[seg_idx]).clamp(0.0, 1.0);
        let ex = x - (ax + t * dx);
        let ey = y - (ay + t * dy);
        ex * ex + ey * ey
    }

    /// Smallest squared distance from (x, y) to segments `start..end`.
    pub fn batch_point_segment_distance_sq(&self, x: f64, y: f64, start: usize, end: usize) -> f64 {
        let mut best = f64::INFINITY;
        for i in start..end {
            let d = self.point_segment_distance_sq(i, x, y);
            if d < best {
                best = d;
            }
        }
        best
    }
}

/// Every obstacle primitive, flattened once per solve.
#[derive(Debug, Clone, Default)]
pub struct ObstacleSet {
    /// Point obstacles in input order of first appearance, exact duplicates
    /// dropped: a repeated site adds no distance information and would only
    /// produce degenerate Voronoi faces.
    pub points: Vec<[f64; 2]>,
    /// Segment obstacles covering every line and polygon edge.
    pub segments: SegmentIndex,
    /// True when `segments` covers a closed area, so the region between its
    /// edges is unreachable without crossing an obstacle.
    pub has_areal_obstacle: bool,
}

impl ObstacleSet {
    /// Flatten any obstacle geometry into its point and segment primitives.
    ///
    /// A geometry collection is flattened recursively; other variants map to
    /// points or segments directly. Empty input is rejected, because there is
    /// no largest empty circle among no obstacles.
    pub fn from_geometry(geom: &Geometry<f64>) -> Result<Self, MecError> {
        let mut set = Self::default();
        set.push_geometry(geom, 0)?;
        if set.points.is_empty() && set.segments.is_empty() {
            return Err(MecError::InvalidInput(
                "obstacle geometry has no points or segments".to_string(),
            ));
        }
        set.validate_finite()?;
        Ok(set)
    }

    fn push_geometry(&mut self, geom: &Geometry<f64>, owner: usize) -> Result<(), MecError> {
        match geom {
            Geometry::Point(p) => self.push_point(p),
            Geometry::MultiPoint(mp) => {
                for p in &mp.0 {
                    self.push_point(p);
                }
            }
            Geometry::Line(l) => self.push_line(l, owner),
            Geometry::LineString(ls) => self.push_line_string(ls, owner)?,
            Geometry::MultiLineString(mls) => {
                for ls in &mls.0 {
                    self.push_line_string(ls, owner)?;
                }
            }
            Geometry::Polygon(p) => self.push_polygon(p, owner)?,
            Geometry::MultiPolygon(mp) => {
                for p in &mp.0 {
                    self.push_polygon(p, owner)?;
                }
            }
            Geometry::Rect(r) => {
                let poly = r.to_polygon();
                self.push_line_string(poly.exterior(), owner)?;
            }
            Geometry::Triangle(t) => {
                let poly = t.to_polygon();
                self.push_line_string(poly.exterior(), owner)?;
            }
            Geometry::GeometryCollection(gc) => {
                for g in &gc.0 {
                    self.push_geometry(g, owner)?;
                }
            }
        }
        Ok(())
    }

    fn push_point(&mut self, p: &Point<f64>) {
        let xy = [p.x(), p.y()];
        if !self.points.contains(&xy) {
            self.points.push(xy);
        }
    }

    fn push_line(&mut self, l: &Line<f64>, owner: usize) {
        self.segments
            .push([l.start.x, l.start.y], [l.end.x, l.end.y], owner);
    }

    fn push_line_string(&mut self, ls: &LineString<f64>, owner: usize) -> Result<(), MecError> {
        if ls.0.len() < 2 {
            return Err(MecError::InvalidInput(
                "obstacle line string has fewer than 2 coordinates".to_string(),
            ));
        }
        for w in ls.0.windows(2) {
            self.segments
                .push([w[0].x, w[0].y], [w[1].x, w[1].y], owner);
        }
        Ok(())
    }

    fn push_polygon(&mut self, poly: &Polygon<f64>, owner: usize) -> Result<(), MecError> {
        self.has_areal_obstacle = true;
        self.push_line_string(poly.exterior(), owner)?;
        for hole in poly.interiors() {
            self.push_line_string(hole, owner)?;
        }
        Ok(())
    }

    /// Reject any non-finite coordinate in either table.
    pub fn validate_finite(&self) -> Result<(), MecError> {
        for p in &self.points {
            if !p[0].is_finite() || !p[1].is_finite() {
                return Err(MecError::InvalidInput(
                    "obstacle point has a NaN or infinite coordinate".to_string(),
                ));
            }
        }
        let s = &self.segments;
        for i in 0..s.len() {
            if !s.ax[i].is_finite()
                || !s.ay[i].is_finite()
                || !s.bx[i].is_finite()
                || !s.by[i].is_finite()
            {
                return Err(MecError::InvalidInput(
                    "obstacle segment has a NaN or infinite coordinate".to_string(),
                ));
            }
        }
        Ok(())
    }

    /// Bounds of every obstacle as (min_x, min_y, max_x, max_y), or None.
    pub fn bounds(&self) -> Option<(f64, f64, f64, f64)> {
        let mut min_x = f64::INFINITY;
        let mut min_y = f64::INFINITY;
        let mut max_x = f64::NEG_INFINITY;
        let mut max_y = f64::NEG_INFINITY;
        let mut any = false;

        for p in &self.points {
            any = true;
            min_x = min_x.min(p[0]);
            min_y = min_y.min(p[1]);
            max_x = max_x.max(p[0]);
            max_y = max_y.max(p[1]);
        }
        for i in 0..self.segments.len() {
            any = true;
            min_x = min_x.min(self.segments.bbox_minx[i]);
            min_y = min_y.min(self.segments.bbox_miny[i]);
            max_x = max_x.max(self.segments.bbox_maxx[i]);
            max_y = max_y.max(self.segments.bbox_maxy[i]);
        }

        any.then_some((min_x, min_y, max_x, max_y))
    }

    /// Smallest squared distance from (x, y) to any obstacle primitive.
    #[inline]
    pub fn clearance_sq(&self, x: f64, y: f64) -> f64 {
        let mut best = f64::INFINITY;
        for p in &self.points {
            let dx = x - p[0];
            let dy = y - p[1];
            let d = dx * dx + dy * dy;
            if d < best {
                best = d;
            }
        }
        for i in 0..self.segments.len() {
            let d = self.segments.point_segment_distance_sq(i, x, y);
            if d < best {
                best = d;
            }
        }
        best
    }

    /// Number of point obstacles.
    pub fn point_count(&self) -> usize {
        self.points.len()
    }

    /// True when every obstacle is a point, so the exact path applies.
    pub fn is_point_only(&self) -> bool {
        !self.points.is_empty() && self.segments.is_empty()
    }

    /// Every obstacle vertex, for the convex-hull default region.
    pub fn vertices(&self) -> Vec<[f64; 2]> {
        let mut out: Vec<[f64; 2]> = self.points.clone();
        for i in 0..self.segments.len() {
            out.push([self.segments.ax[i], self.segments.ay[i]]);
            out.push([self.segments.bx[i], self.segments.by[i]]);
        }
        out
    }
}

/// The region the circle centre must lie in, if the caller supplied one.
#[derive(Debug, Clone, Default)]
pub struct BoundarySet {
    /// Boundary rings as closed coordinate loops, one entry per ring.
    pub rings: Vec<Vec<[f64; 2]>>,
}

impl BoundarySet {
    /// Flatten a polygonal boundary into rings.
    ///
    /// Only polygonal input is accepted, because the boundary constrains
    /// containment and "inside" is undefined for a point or a line.
    pub fn from_polygonal(geom: &Geometry<f64>) -> Result<Self, MecError> {
        let mut set = Self::default();
        match geom {
            Geometry::Polygon(p) => set.push_polygon(p),
            Geometry::MultiPolygon(MultiPolygon(mp)) => {
                for p in mp {
                    set.push_polygon(p);
                }
            }
            Geometry::Rect(r) => {
                let poly = r.to_polygon();
                set.push_polygon(&poly);
            }
            Geometry::Triangle(t) => {
                let poly = t.to_polygon();
                set.push_polygon(&poly);
            }
            other => {
                return Err(MecError::InvalidInput(format!(
                    "boundary must be polygonal, got {}",
                    variant_name(other)
                )));
            }
        }
        if set.rings.is_empty() {
            return Err(MecError::InvalidInput("boundary is empty".to_string()));
        }
        Ok(set)
    }

    fn push_polygon(&mut self, poly: &Polygon<f64>) {
        let mut push_ring = |ls: &LineString<f64>| {
            let ring: Vec<[f64; 2]> = ls.0.iter().map(|c| [c.x, c.y]).collect();
            if ring.len() >= 4 {
                self.rings.push(ring);
            }
        };
        push_ring(poly.exterior());
        for hole in poly.interiors() {
            push_ring(hole);
        }
    }

    /// Bounds of the boundary as (min_x, min_y, max_x, max_y), or None.
    pub fn bounds(&self) -> Option<(f64, f64, f64, f64)> {
        let mut min_x = f64::INFINITY;
        let mut min_y = f64::INFINITY;
        let mut max_x = f64::NEG_INFINITY;
        let mut max_y = f64::NEG_INFINITY;
        let mut any = false;
        for ring in &self.rings {
            for p in ring {
                any = true;
                min_x = min_x.min(p[0]);
                min_y = min_y.min(p[1]);
                max_x = max_x.max(p[0]);
                max_y = max_y.max(p[1]);
            }
        }
        any.then_some((min_x, min_y, max_x, max_y))
    }
}

/// The `Geometry` variant name, for error messages.
pub(crate) fn variant_name(geom: &Geometry<f64>) -> &'static str {
    match geom {
        Geometry::Point(_) => "Point",
        Geometry::Line(_) => "Line",
        Geometry::LineString(_) => "LineString",
        Geometry::Polygon(_) => "Polygon",
        Geometry::MultiPoint(_) => "MultiPoint",
        Geometry::MultiLineString(_) => "MultiLineString",
        Geometry::MultiPolygon(_) => "MultiPolygon",
        Geometry::GeometryCollection(_) => "GeometryCollection",
        Geometry::Rect(_) => "Rect",
        Geometry::Triangle(_) => "Triangle",
    }
}

/// Build a closed line string from an open pair list.
pub fn coords_to_line_string(pts: &[(f64, f64)]) -> LineString<f64> {
    LineString::from(
        pts.iter()
            .map(|p| Coord { x: p.0, y: p.1 })
            .collect::<Vec<_>>(),
    )
}
