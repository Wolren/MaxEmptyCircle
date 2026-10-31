//! Exact refinement of the grid path's answer.
//!
//! The grid ranks cells by a bound, so its incumbent is the best *sampled*
//! centre, and it can sit a fraction short of the true peak. Walking the field
//! cannot close that reliably: near the peak the clearance is flat to first order,
//! so the gradient carries almost no signal, and at a boundary optimum the only
//! improving direction runs along the region edge, which no compass step can take.
//!
//! So this module does not walk. It enumerates, in two places: the interior
//! three-way apex here, and the boundary maximum per region edge below.
//!
//! # Why a triple enumeration is complete
//!
//! Let `f(p) = min_i d_i(p)` be the clearance and suppose `p` is an interior
//! local maximum. Each `d_i` is smooth at `p` once the nearest point of
//! primitive `i` is fixed, with `grad d_i` a unit vector: `(p - q_i)/|p - q_i|`
//! for a point `q_i`, and `+/- n_i` for a segment whose nearest point is
//! interior, `n_i` the unit normal of its supporting line. A local maximum
//! needs `0` in the subdifferential, which is the convex hull of those unit
//! vectors. In the plane Carathéodory's theorem then gives a subset of at most
//! three of them whose hull already contains `0`, so at least three distances
//! are active and equal, or two are active with opposite gradients.
//!
//! The two-vector case is `grad d_1 = -grad d_2`. For two points that puts `p`
//! between them on the line through them, a measure-zero coincidence; for two
//! segments it needs parallel supporting lines, and there the clearance is
//! exactly constant along the whole midline, so every point of it is optimal and
//! no refinement is needed. Every other interior maximum is a point equidistant
//! from three primitives whose nearest points are three of the primitive's own
//! *carriers*: a point carrier, or the supporting line of a segment whose
//! nearest point falls inside it.
//!
//! So: enumerate the carriers near the incumbent, solve the equidistance system
//! for each triple in closed form, and keep the best point that is admissible.
//! That set contains the true apex, and it is solved exactly rather than
//! approached.
//!
//! # Why the carriers near the incumbent are enough
//!
//! All three active primitives are within `f(p)` of the apex, and the incumbent
//! is within one grid tolerance of it, so they are all within a small relative
//! band of the incumbent's clearance. Taking the band wider than the grid
//! tolerance and capping the count keeps the cost at a few hundred triples while
//! containing every triple that can form the apex.
//!
//! A candidate that does not win is simply dropped, so this stage can never make
//! an answer worse than the incumbent it was given.
use super::super::boundary::Region;
use super::super::index::NearestObstacleIndex;
use super::super::input::ObstacleSet;
use super::exact::circumcentre;
use crate::tuning;

/// A convex carrier whose distance function is locally the primitive's distance.
#[derive(Debug, Clone, Copy)]
enum Carrier {
    /// A point obstacle or a segment endpoint, from primitive `prim`.
    Point { x: f64, y: f64, prim: u32 },
    /// The supporting line of segment `seg`, as `nx*x + ny*y - c` with
    /// `(nx, ny)` a unit normal. Its distance is the absolute value, and it only
    /// equals the segment's distance while the perpendicular foot stays inside
    /// the segment.
    Line {
        nx: f64,
        ny: f64,
        c: f64,
        seg: u32,
        prim: u32,
    },
}

impl Carrier {
    /// Ordinal of the primitive this carrier belongs to.
    #[inline]
    fn prim(&self) -> u32 {
        match *self {
            Carrier::Point { prim, .. } | Carrier::Line { prim, .. } => prim,
        }
    }
}

/// A point and a direction (not necessarily unit): the line `p0 + t * dir`.
type ParamLine = ([f64; 2], [f64; 2]);

/// One region edge: the base point and the direction its parameter runs along.
#[derive(Debug, Clone, Copy)]
struct EdgeGeom {
    ax: f64,
    ay: f64,
    ex: f64,
    ey: f64,
}

impl EdgeGeom {
    /// Point at edge parameter `t`.
    #[inline]
    fn point(&self, t: f64) -> (f64, f64) {
        (self.ax + t * self.ex, self.ay + t * self.ey)
    }

    /// Edge length.
    #[inline]
    fn len(&self) -> f64 {
        self.ex.hypot(self.ey)
    }

    /// Squared edge length.
    #[inline]
    fn len_sq(&self) -> f64 {
        self.ex * self.ex + self.ey * self.ey
    }
}

/// Refine the incumbent towards the exact apex of its own basin.
///
/// `start` is `(x, y, distance)`. The result is never worse than `start`, and is
/// the exact three-way apex when one exists within the carrier band.
pub fn refine_apex(
    set: &ObstacleSet,
    nb: &NearestObstacleIndex,
    region: &Region,
    start: (f64, f64, f64),
) -> (f64, f64, f64) {
    if !start.0.is_finite() || !start.1.is_finite() || !start.2.is_finite() || start.2 <= 0.0 {
        return start;
    }

    // Primitive ordinals within a small relative band of the incumbent, nearest
    // first so the cap keeps the ones that can support the apex.
    let band = 1.0 + tuning::MEC_APEX_BAND_REL;
    let reach_sq = start.2 * start.2 * band * band;
    let prims = nb.nearest_primitives(start.0, start.1, reach_sq, tuning::MEC_APEX_MAX_PRIMS);
    if prims.len() < 3 {
        return start;
    }

    let carriers = carriers_of(set, &prims);
    let n = carriers.len();
    if n < 3 {
        return start;
    }

    let mut best = start;
    let ctx = ApexCtx {
        set,
        nb,
        region,
        prims: &prims,
    };
    let mut candidates: Vec<[f64; 2]> = Vec::with_capacity(8);
    for i in 0..n {
        for j in (i + 1)..n {
            // Two carriers of one primitive cannot both be active: a primitive's
            // nearest point is unique, so one of the two is not achieved.
            if carriers[i].prim() == carriers[j].prim() {
                continue;
            }
            for k in (j + 1)..n {
                if carriers[k].prim() == carriers[i].prim()
                    || carriers[k].prim() == carriers[j].prim()
                {
                    continue;
                }
                candidates.clear();
                let triple = [carriers[i], carriers[j], carriers[k]];
                solve_three(&triple[0], &triple[1], &triple[2], &mut candidates);
                for c in &candidates {
                    consider(&ctx, &mut best, &triple, c[0], c[1]);
                }
            }
        }
    }
    best
}

/// Refine the answer along the region boundary.
///
/// A constrained optimum sits on a region edge whenever the clearance is still
/// increasing outward there, and on that edge only one degree of freedom is left.
/// A 2-D walk cannot follow it: the only improving direction is *exactly* along
/// the edge, every direction a compass can step in either leaves the region or
/// moves inward, and the walk stops a fraction short. Two obstacles tied along one
/// edge — the common case, since a gap between two obstacles ends where the region
/// does — then costs a relative per-mille that ascent never recovers.
///
/// So the boundary is solved rather than searched.
///
/// # Why the bisector crossings are the whole answer
///
/// Restricted to one edge, the clearance is `g(t) = min_i d_i(t)`, and every
/// `d_i` is the distance to a convex primitive along a line, which is convex in
/// `t`. Between two consecutive points where the nearest primitive changes, `g`
/// equals a single convex `d_i`, and a convex function takes its maximum on a
/// closed interval at an endpoint. So the maximum of the edge is at a point where
/// two primitives are equidistant, or at an edge end. Enumerating the crossings
/// closes the family exactly, in closed form, with no assumption about the shape
/// between samples.
///
/// # Why sampling is still there
///
/// Two jobs, both cheap. It supplies the incumbent so a crossing only has to be
/// solved where it could win, and it selects which primitives can be active: the
/// peak is within half a sample gap of its nearest sample by the 1-Lipschitz
/// property, so a primitive active at the peak is within two gaps of that sample's
/// clearance. Everything farther cannot be one of the two.
pub fn refine_boundary(
    set: &ObstacleSet,
    nb: &NearestObstacleIndex,
    region: &Region,
    start: (f64, f64, f64),
) -> (f64, f64, f64) {
    if !start.0.is_finite() || !start.1.is_finite() || !start.2.is_finite() || start.2 <= 0.0 {
        return start;
    }
    let steps = tuning::MEC_EDGE_SAMPLES;
    if steps == 0 || region.ring_count() == 0 {
        return start;
    }

    /// One edge with the bound its samples put on the whole edge.
    struct Edge {
        geom: EdgeGeom,
        bound: f64,
    }

    let mut best = start;
    let mut edges: Vec<Edge> = Vec::new();

    for r in 0..region.ring_count() {
        let ring = region.ring_coords(r);
        for w in ring.windows(2) {
            let geom = EdgeGeom {
                ax: w[0][0],
                ay: w[0][1],
                ex: w[1][0] - w[0][0],
                ey: w[1][1] - w[0][1],
            };
            let len = geom.len();
            if !len.is_finite() || len <= 0.0 {
                continue;
            }
            let mut sampled = f64::NEG_INFINITY;
            for k in 0..=steps {
                let (x, y) = geom.point(k as f64 / steps as f64);
                let d = clearance(nb, x, y);
                if d > sampled {
                    sampled = d;
                    if d > best.2 {
                        best = (x, y, d);
                    }
                }
            }
            // No point of the edge can clear more than the sampled maximum plus
            // half a sample gap, because the clearance is 1-Lipschitz.
            edges.push(Edge {
                geom,
                bound: sampled + 0.5 * len / steps as f64,
            });
        }
    }

    // Bound order, not sampled order: a narrow peak loses up to half a gap
    // between samples, so the edge holding the true maximum can sample below an
    // edge that merely samples well. The bound also gives the stopping rule a
    // proof rather than a guess.
    edges.sort_by(|a, b| {
        b.bound
            .partial_cmp(&a.bound)
            .unwrap_or(std::cmp::Ordering::Equal)
    });
    for edge in edges.iter() {
        if edge.bound <= best.2 {
            break;
        }
        best = refine_edge(set, nb, &edge.geom, steps, best);
    }
    best
}

/// Exact maximum of the clearance along one edge.
///
/// `start` is the incumbent; the result is never worse.
fn refine_edge(
    set: &ObstacleSet,
    nb: &NearestObstacleIndex,
    edge: &EdgeGeom,
    steps: usize,
    start: (f64, f64, f64),
) -> (f64, f64, f64) {
    let len = edge.len();
    if !len.is_finite() || len <= 0.0 || steps == 0 {
        return start;
    }
    let e2 = edge.len_sq();
    let half_gap = 0.5 * len / steps as f64;

    let mut best = start;
    let mut values: Vec<f64> = Vec::with_capacity(steps + 1);
    for k in 0..=steps {
        let (x, y) = edge.point(k as f64 / steps as f64);
        let d = clearance(nb, x, y);
        values.push(d);
        if d > best.2 {
            best = (x, y, d);
        }
    }

    // Primitives that can be active at a peak on this edge. A sample whose own
    // clearance plus a gap already loses cannot hide one, and the active distance
    // at the peak is at most the best sample plus two gaps.
    let mut ords: Vec<usize> = Vec::new();
    for (k, &v) in values.iter().enumerate() {
        if v + half_gap <= best.2 {
            continue;
        }
        let (x, y) = edge.point(k as f64 / steps as f64);
        let reach = v + 3.0 * half_gap;
        ords.extend(nb.supporting_primitives(x, y, reach * reach));
    }
    if ords.len() < 2 {
        return best;
    }
    ords.sort_unstable();
    ords.dedup();

    let cells: Vec<(f64, usize)> = ords.iter().map(|&i| (0.0, i)).collect();
    let carriers = carriers_of(set, &cells);
    if carriers.len() < 2 {
        return best;
    }
    let edge_carriers: Vec<EdgeCarrier> =
        carriers.iter().map(|c| EdgeCarrier::new(c, edge)).collect();

    let mut roots: Vec<f64> = Vec::with_capacity(2);
    for i in 0..edge_carriers.len() {
        for j in (i + 1)..edge_carriers.len() {
            roots.clear();
            edge_carriers[i].equal_roots(&edge_carriers[j], e2, &mut roots);
            for &t in &roots {
                if !(-1e-9..=1.0 + 1e-9).contains(&t) {
                    continue;
                }
                let t = t.clamp(0.0, 1.0);
                // The nominal radius is the common distance at the crossing. It is
                // a cheap upper bound on the measured clearance, so a crossing that
                // cannot beat the incumbent never pays for a hierarchy query.
                let nominal = edge_carriers[i].distance_sq(t, e2);
                if nominal <= best.2 * best.2 {
                    continue;
                }
                let (x, y) = edge.point(t);
                let d = clearance(nb, x, y);
                if d > best.2 {
                    best = (x, y, d);
                }
            }
        }
    }
    best
}
/// A carrier's squared distance along an edge, as a function of the edge parameter.
#[derive(Debug, Clone, Copy)]
enum EdgeCarrier {
    /// A point: `|E|^2 t^2 + b t + c`, with `|E|^2` shared by every point carrier
    /// on the edge and therefore passed in rather than stored.
    Point { b: f64, c: f64 },
    /// A line: `(alpha t + beta)^2`.
    Line { alpha: f64, beta: f64 },
}

impl EdgeCarrier {
    fn new(c: &Carrier, edge: &EdgeGeom) -> Self {
        match *c {
            Carrier::Point { x, y, .. } => {
                let dx = edge.ax - x;
                let dy = edge.ay - y;
                Self::Point {
                    b: 2.0 * (edge.ex * dx + edge.ey * dy),
                    c: dx * dx + dy * dy,
                }
            }
            Carrier::Line { nx, ny, c, .. } => Self::Line {
                alpha: nx * edge.ex + ny * edge.ey,
                beta: nx * edge.ax + ny * edge.ay - c,
            },
        }
    }

    /// Squared distance at edge parameter `t`.
    #[inline]
    fn distance_sq(&self, t: f64, e2: f64) -> f64 {
        match *self {
            Self::Point { b, c } => (e2 * t + b) * t + c,
            Self::Line { alpha, beta } => {
                let s = alpha * t + beta;
                s * s
            }
        }
    }

    /// Parameters where this carrier and `other` are equidistant.
    fn equal_roots(&self, other: &Self, e2: f64, out: &mut Vec<f64>) {
        match (*self, *other) {
            // The quadratic terms are equal, so the difference is linear.
            (Self::Point { b: b1, c: c1 }, Self::Point { b: b2, c: c2 }) => {
                let db = b1 - b2;
                let dc = c1 - c2;
                if db != 0.0 && db.is_finite() {
                    out.push(-dc / db);
                }
            }
            (Self::Point { b, c }, Self::Line { alpha, beta }) => {
                // (e2 - alpha^2) t^2 + (b - 2 alpha beta) t + (c - beta^2) = 0
                quadratic(
                    e2 - alpha * alpha,
                    b - 2.0 * alpha * beta,
                    c - beta * beta,
                    e2,
                    out,
                );
            }
            (Self::Line { alpha, beta }, Self::Point { b, c }) => {
                quadratic(
                    e2 - alpha * alpha,
                    b - 2.0 * alpha * beta,
                    c - beta * beta,
                    e2,
                    out,
                );
            }
            // A difference of squares: the two angle bisectors, one per sign.
            (
                Self::Line {
                    alpha: a1,
                    beta: b1,
                },
                Self::Line {
                    alpha: a2,
                    beta: b2,
                },
            ) => {
                let d = a1 - a2;
                if d != 0.0 && d.is_finite() {
                    out.push(-(b1 - b2) / d);
                }
                let s = a1 + a2;
                if s != 0.0 && s.is_finite() {
                    out.push(-(b1 + b2) / s);
                }
            }
        }
    }
}

/// Real roots of `a t^2 + b t + c`, with the leading coefficient judged
/// degenerate against the edge's own scale.
fn quadratic(a: f64, b: f64, c: f64, e2: f64, out: &mut Vec<f64>) {
    if a.abs() <= 1e-14 * e2 {
        if b != 0.0 && b.is_finite() {
            out.push(-c / b);
        }
        return;
    }
    let mut disc = b * b - 4.0 * a * c;
    if disc < 0.0 {
        if disc < -1e-12 * (b * b + (4.0 * a * c).abs()) {
            return;
        }
        disc = 0.0;
    }
    let root = disc.sqrt();
    out.push((-b - root) / (2.0 * a));
    let t2 = (-b + root) / (2.0 * a);
    if (t2 - out[out.len() - 1]).abs() > 1e-12 * (1.0 + t2.abs()) {
        out.push(t2);
    }
}

/// Clearance at a point, floored at zero.
#[inline]
fn clearance(nb: &NearestObstacleIndex, x: f64, y: f64) -> f64 {
    match nb.nearest_distance_sq(x, y) {
        Some((d, _)) if d.is_finite() && d > 0.0 => d.sqrt(),
        _ => 0.0,
    }
}

/// The carriers of the given primitive ordinals.
///
/// A segment contributes three: its supporting line, which covers a nearest point
/// inside the segment, and both endpoints, which cover a nearest point at either
/// end. Point obstacles contribute themselves.
fn carriers_of(set: &ObstacleSet, prims: &[(f64, usize)]) -> Vec<Carrier> {
    let n_points = set.points.len();
    let mut out: Vec<Carrier> = Vec::with_capacity(prims.len() * 3);
    for &(_, prim) in prims {
        if prim < n_points {
            let p = set.points[prim];
            out.push(Carrier::Point {
                x: p[0],
                y: p[1],
                prim: prim as u32,
            });
            continue;
        }
        let s = prim - n_points;
        if s >= set.segments.len() {
            continue;
        }
        let (ax, ay) = (set.segments.ax[s], set.segments.ay[s]);
        let (bx, by) = (set.segments.bx[s], set.segments.by[s]);
        let (dx, dy) = (bx - ax, by - ay);
        let len = dx.hypot(dy);
        if len > 0.0 {
            let (nx, ny) = (-dy / len, dx / len);
            out.push(Carrier::Line {
                nx,
                ny,
                c: nx * ax + ny * ay,
                seg: s as u32,
                prim: prim as u32,
            });
        }
        out.push(Carrier::Point {
            x: ax,
            y: ay,
            prim: prim as u32,
        });
        out.push(Carrier::Point {
            x: bx,
            y: by,
            prim: prim as u32,
        });
    }
    out
}

/// Keep an admissible candidate that beats the incumbent.
///
/// Two filters run before the region test and the hierarchy query, and both are
/// sound because they only reject points that cannot be the maximum.
///
/// The first is the carrier test: every primitive's distance is at least its
/// carrier's distance, so a candidate the incumbent already clears is not worth
/// measuring. The primitives are tested rather than the carriers, because a
/// carrier is *not* an upper bound on its primitive: a segment's supporting line
/// is nearer than the segment whenever the perpendicular foot falls outside it.
///
/// The second is carrier validity. A maximum is equidistant from three
/// primitives through three *achieved* carriers: a segment's supporting line
/// counts only where its foot lands inside the segment, and otherwise the active
/// carrier is an endpoint, which is enumerated separately. A triple whose foot
/// falls outside describes a point where the nominal radius is not reached at
/// all, so it cannot be the maximum and is dropped. This is what keeps a box's
/// four bounding lines from producing a candidate per corner pair.
/// What a candidate is scored against.
struct ApexCtx<'a> {
    set: &'a ObstacleSet,
    nb: &'a NearestObstacleIndex,
    region: &'a Region,
    prims: &'a [(f64, usize)],
}

#[inline]
fn consider(ctx: &ApexCtx, best: &mut (f64, f64, f64), triple: &[Carrier; 3], x: f64, y: f64) {
    if !x.is_finite() || !y.is_finite() {
        return;
    }
    let bar = best.2 * best.2;
    for &(_, prim) in ctx.prims {
        if ctx.nb.prim_distance_sq(prim, x, y) < bar {
            return;
        }
    }
    for car in triple {
        if let Carrier::Line { seg, .. } = *car {
            if !foot_inside(ctx.set, seg as usize, x, y) {
                return;
            }
        }
    }
    if !ctx.region.contains_closed(x, y) {
        return;
    }
    let Some((d_sq, _)) = ctx.nb.nearest_distance_sq(x, y) else {
        return;
    };
    if d_sq > bar {
        *best = (x, y, d_sq.sqrt());
    }
}

/// True when the perpendicular foot of (x, y) on segment `s` lies inside it.
#[inline]
fn foot_inside(set: &ObstacleSet, s: usize, x: f64, y: f64) -> bool {
    let seg = &set.segments;
    let t = ((x - seg.ax[s]) * seg.dir_x[s] + (y - seg.ay[s]) * seg.dir_y[s]) * seg.inv_len_sq[s];
    // A hair of slack: a foot a rounding error outside the segment is the
    // endpoint case, which the endpoint carrier already covers.
    (-1e-9..=1.0 + 1e-9).contains(&t)
}

/// Candidate points equidistant from three carriers.
fn solve_three(a: &Carrier, b: &Carrier, c: &Carrier, out: &mut Vec<[f64; 2]>) {
    let mut pts = [[0.0f64; 2]; 3];
    let mut lines = [[0.0f64; 3]; 3];
    let mut n_pts = 0usize;
    let mut n_lines = 0usize;
    for car in [a, b, c] {
        match *car {
            Carrier::Point { x, y, .. } => {
                pts[n_pts] = [x, y];
                n_pts += 1;
            }
            Carrier::Line { nx, ny, c, .. } => {
                lines[n_lines] = [nx, ny, c];
                n_lines += 1;
            }
        }
    }

    match n_lines {
        // Three points: their circumcentre.
        0 => {
            if let Some((x, y)) = circumcentre(pts[0], pts[1], pts[2]) {
                out.push([x, y]);
            }
        }
        // Two points and a line: the centre sits on the points' perpendicular
        // bisector, and the remaining condition is quadratic along it.
        1 => {
            if let Some((p0, dir)) = perp_bisector(pts[0], pts[1]) {
                point_line_roots(p0, dir, lines[0], pts[0], out);
            }
        }
        // A point and two lines: each sign pair of the lines gives one angle
        // bisector, and the condition against the point is quadratic along it.
        2 => {
            for s0 in [1.0, -1.0] {
                for s1 in [1.0, -1.0] {
                    let Some((p0, dir)) = angle_bisector(lines[0], s0, lines[1], s1) else {
                        continue;
                    };
                    let signed = [s0 * lines[0][0], s0 * lines[0][1], s0 * lines[0][2]];
                    point_line_roots(p0, dir, signed, pts[0], out);
                }
            }
        }
        // Three lines: intersect two angle bisectors, one per pair. Fixing the
        // first line's sign covers every solution, because flipping all three
        // signs leaves each bisector where it was.
        _ => {
            for s1 in [1.0, -1.0] {
                for s2 in [1.0, -1.0] {
                    let Some(first) = angle_bisector(lines[0], 1.0, lines[1], s1) else {
                        continue;
                    };
                    let Some(second) = angle_bisector(lines[0], 1.0, lines[2], s2) else {
                        continue;
                    };
                    if let Some(p) = intersect(first, second) {
                        out.push(p);
                    }
                }
            }
        }
    }
}

/// Perpendicular bisector of two points, as a parameterised line.
///
/// The direction is not normalised, and nothing below needs it to be: the
/// quadratic form carries the direction's length explicitly.
fn perp_bisector(a: [f64; 2], b: [f64; 2]) -> Option<ParamLine> {
    let (dx, dy) = (b[0] - a[0], b[1] - a[1]);
    let len_sq = dx * dx + dy * dy;
    if !len_sq.is_finite() || len_sq <= 0.0 {
        return None;
    }
    let mid = [(a[0] + b[0]) * 0.5, (a[1] + b[1]) * 0.5];
    Some((mid, [-dy, dx]))
}

/// Bisector of two lines with orientation signs, as a parameterised line.
///
/// The bisector is `{p : s0*(n0.p - c0) = s1*(n1.p - c1)}`, a line with normal
/// `N = s0*n0 - s1*n1` and offset `C = s0*c0 - s1*c1`, so `N.p = C`. The
/// direction is left unnormalised; the caller's quadratic form does not need a
/// unit vector, and skipping the square root keeps a triple solve at a handful
/// of multiplies.
fn angle_bisector(l0: [f64; 3], s0: f64, l1: [f64; 3], s1: f64) -> Option<ParamLine> {
    let nx = s0 * l0[0] - s1 * l1[0];
    let ny = s0 * l0[1] - s1 * l1[1];
    let c = s0 * l0[2] - s1 * l1[2];
    let len_sq = nx * nx + ny * ny;
    // Parallel or coincident lines: no bisector, or every point of it.
    if !len_sq.is_finite() || len_sq <= 1e-24 {
        return None;
    }
    let p0 = [nx * c / len_sq, ny * c / len_sq];
    Some((p0, [-ny, nx]))
}

/// Points on `p0 + t*dir` that are equidistant from a point and a line.
///
/// `signed` is the line with the sign that makes `n.p - c` positive on the side
/// the solution lies; both signs produce the same squared equation, so the
/// squaring below loses nothing.
fn point_line_roots(
    p0: [f64; 2],
    dir: [f64; 2],
    signed: [f64; 3],
    q: [f64; 2],
    out: &mut Vec<[f64; 2]>,
) {
    let (nx, ny, c) = (signed[0], signed[1], signed[2]);
    let (dx, dy) = (dir[0], dir[1]);
    let nd = nx * dx + ny * dy;
    let dd = dx * dx + dy * dy;

    let rx = p0[0] - q[0];
    let ry = p0[1] - q[1];
    let s0 = nx * p0[0] + ny * p0[1] - c;

    // |p - q|^2 = (n.p - c)^2 with p = p0 + t*dir, dir not normalised:
    //   a t^2 + b t + c0 = 0
    //   a  = |dir|^2 - (n.dir)^2   (>0 unless dir is along the line's normal)
    //   b  = 2 * (dir.(p0 - q) - s0 * n.dir)
    //   c0 = |p0 - q|^2 - s0^2
    let a = dd - nd * nd;
    let b = 2.0 * (dx * rx + dy * ry - s0 * nd);
    let c0 = rx * rx + ry * ry - s0 * s0;

    // a scales with |dir|^2, so the degeneracy test is relative: a is a fraction
    // |dir|^2*sin^2(angle) of the direction's own length.
    if a <= dd * 1e-24 {
        if b.abs() <= 1e-12 * dd {
            return;
        }
        push_on_line(p0, dir, -c0 / b, out);
        return;
    }
    let mut disc = b * b - 4.0 * a * c0;
    if disc < 0.0 {
        // A double root landing a rounding error below zero is still a root.
        let scale = b * b + (4.0 * a * c0).abs();
        if disc < -1e-12 * scale {
            return;
        }
        disc = 0.0;
    }
    let root = disc.sqrt();
    let t1 = (-b - root) / (2.0 * a);
    let t2 = (-b + root) / (2.0 * a);
    push_on_line(p0, dir, t1, out);
    if (t2 - t1).abs() > 1e-12 * (1.0 + t1.abs()) {
        push_on_line(p0, dir, t2, out);
    }
}

#[inline]
fn push_on_line(p0: [f64; 2], dir: [f64; 2], t: f64, out: &mut Vec<[f64; 2]>) {
    let p = [p0[0] + t * dir[0], p0[1] + t * dir[1]];
    if p[0].is_finite() && p[1].is_finite() {
        out.push(p);
    }
}

/// Intersection of two parameterised lines.
///
/// The directions are not normalised, so the parallel test is relative to their
/// lengths: `|d x e|` is `|d||e|sin(angle)`.
fn intersect(a: ParamLine, b: ParamLine) -> Option<[f64; 2]> {
    let (p, d) = a;
    let (q, e) = b;
    let den = d[0] * e[1] - d[1] * e[0];
    let scale = (d[0] * d[0] + d[1] * d[1]) * (e[0] * e[0] + e[1] * e[1]);
    if !den.is_finite() || den * den <= 1e-24 * scale {
        return None;
    }
    let t = ((q[0] - p[0]) * e[1] - (q[1] - p[1]) * e[0]) / den;
    let x = p[0] + t * d[0];
    let y = p[1] + t * d[1];
    (x.is_finite() && y.is_finite()).then_some([x, y])
}
