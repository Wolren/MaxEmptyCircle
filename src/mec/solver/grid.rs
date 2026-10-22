//! Grid solver: tiles the search box with cells and refines the best one.
//!
//! The search is a branch and bound on one invariant: the distance to the
//! nearest obstacle is 1-Lipschitz, so for a cell of half-size `h` centred at
//! `c`, no point of that cell can be farther from an obstacle than
//! `obstacle_distance(c) + h * sqrt(2)`. That bound holds for every centre,
//! which is what makes the pruning sound for cells that straddle an obstacle.
//! A cell leaves the queue only when its bound proves it cannot beat the
//! incumbent by more than the tolerance.
//!
//! The empty circle differs from the inscribed circle in one way that shapes
//! this loop. The circle centre need not be inside any polygon: the obstacles
//! are what it must avoid, and the only region constraint is [`Region`]. So the
//! primal asks "is this centre clear of every obstacle", and the region asks
//! "is this centre allowed here" as a second, separate test.
use std::collections::BinaryHeap;

use super::super::boundary::Region;
use super::super::index::NearestObstacleIndex;
use super::super::input::ObstacleSet;
use crate::tuning;

const SQRT2: f64 = std::f64::consts::SQRT_2;

#[derive(Debug, Clone, Copy)]
struct GridCell {
    x: f64,
    y: f64,
    h_size: f64,
    max_dist: f64,
    /// Primitive ordinal achieving the nearest distance at this cell's centre,
    /// used to seed the next level so it starts from a near hit.
    support: usize,
}

impl PartialEq for GridCell {
    fn eq(&self, other: &Self) -> bool {
        self.max_dist == other.max_dist
    }
}

impl Eq for GridCell {}

impl PartialOrd for GridCell {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for GridCell {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        self.max_dist
            .partial_cmp(&other.max_dist)
            .unwrap_or(std::cmp::Ordering::Equal)
    }
}

/// One grid answer, before it is dressed as a public result.
#[derive(Debug, Clone, Copy)]
pub struct GridSolution {
    /// Centre x.
    pub x: f64,
    /// Centre y.
    pub y: f64,
    /// Squared clearance at the centre.
    pub radius_sq: f64,
    /// Cells the search popped.
    pub iterations: usize,
    /// True when the budget, not the bound, ended the search.
    pub exhausted_budget: bool,
}

/// Solve the largest empty circle by branch and bound over the region.
///
/// The radius of the returned result is measured against the obstacle set at
/// the returned centre rather than read off the cell bound.
#[cfg_attr(feature = "hotpath", hotpath::measure)]
pub fn solve_grid(
    set: &ObstacleSet,
    nb: &NearestObstacleIndex,
    region: &Region,
    tolerance: f64,
) -> Option<GridSolution> {
    let (min_x, min_y, max_x, max_y) = region.bounds();

    let width = max_x - min_x;
    let height = max_y - min_y;
    let cell_size = width.max(height);

    if cell_size <= 0.0 || !tolerance.is_finite() || tolerance <= 0.0 {
        // A collapsed region leaves one candidate: the obstacle centroid.
        let (cx, cy) = centroid_of(set)?;
        if region.contains(cx, cy) {
            let r_sq = set.clearance_sq(cx, cy);
            if r_sq.is_finite() && r_sq > 0.0 {
                return Some(GridSolution {
                    x: cx,
                    y: cy,
                    radius_sq: r_sq,
                    iterations: 0,
                    exhausted_budget: false,
                });
            }
        }
        return None;
    }

    // Half-size of a first-level cell: the box is covered by at most GRID_SIDE
    // cells along its longer axis and by however many are needed along the
    // shorter one, so a thin region is still covered end to end.
    let side = tuning::MEC_GRID_SIDE;
    let h_size = cell_size / (2.0 * side as f64);
    let step = 2.0 * h_size;
    let cells_x = ((width / step).ceil() as usize).max(1);
    let cells_y = ((height / step).ceil() as usize).max(1);

    let mut queue: BinaryHeap<GridCell> = BinaryHeap::with_capacity(cells_x * cells_y * 2);
    let mut best_x = 0.0;
    let mut best_y = 0.0;
    let mut best_dist = 0.0;
    let mut found = false;

    for i in 0..cells_x {
        let cx = min_x + h_size + i as f64 * step;
        for j in 0..cells_y {
            let cy = min_y + h_size + j as f64 * step;
            let d = nearest_distance(nb, cx, cy);
            if region.contains(cx, cy) && d > best_dist {
                best_x = cx;
                best_y = cy;
                best_dist = d;
                found = true;
            }
            queue.push(GridCell {
                x: cx,
                y: cy,
                h_size,
                max_dist: d + h_size * SQRT2,
                support: usize::MAX,
            });
        }
    }

    let mut iterations = 0usize;
    // True when the loop stopped on the budget instead of on the bound, which
    // means the incumbent may still be short of the peak.
    let mut exhausted_budget = false;

    while let Some(cell) = queue.pop() {
        iterations += 1;

        if iterations > tuning::MEC_GRID_CELL_BUDGET {
            exhausted_budget = true;
            break;
        }

        // The heap is ordered by the bound, so once the best remaining cell
        // cannot beat the incumbent by more than the prune gap, no cell can.
        let radius_bound = best_dist + tolerance;
        if found && cell.max_dist <= radius_bound {
            break;
        }

        let half = cell.h_size * 0.5;

        // Refining below the tolerance gains nothing: the child bound differs
        // from the parent's by at most `half * sqrt(2)`, already under the gap
        // the search resolves. Without this floor a cell in a region with no
        // admissible centre has no incumbent to prune it and splits until the
        // half-size underflows.
        if half * SQRT2 <= tolerance {
            continue;
        }

        let children = [
            (cell.x - half, cell.y - half),
            (cell.x + half, cell.y - half),
            (cell.x - half, cell.y + half),
            (cell.x + half, cell.y + half),
        ];

        // A child is pruned when `d + half*sqrt(2) <= best + tolerance`, which
        // rearranges to `d <= best - max(0, half*sqrt(2) - tolerance)`. Asking
        // the index for "the nearest distance, or a value under this bar" lets
        // the descent stop at the first primitive inside the bar instead of
        // proving the true minimum.
        let slack = half * SQRT2 - (radius_bound - best_dist);
        let threshold = if found && slack > 0.0 {
            Some(best_dist - slack)
        } else {
            None
        };

        for (cx, cy) in children {
            // A child wholly outside the search box holds no candidate.
            if cx + half < min_x || cx - half > max_x || cy + half < min_y || cy - half > max_y {
                continue;
            }
            let (d, support) = match threshold {
                Some(t) => match nb.nearest_distance_sq_above(
                    cx,
                    cy,
                    (t.max(0.0)) * (t.max(0.0)),
                    Some(cell.support),
                ) {
                    Some(v) => v,
                    // Proven within the bar: this child cannot improve.
                    None => continue,
                },
                None => {
                    let Some(v) = nb.nearest_distance_sq(cx, cy) else {
                        continue;
                    };
                    v
                }
            };
            let d = d.max(0.0).sqrt();
            let max_dist = d + half * SQRT2;
            // The bound test needs no containment result, so it runs first and
            // most cells never pay for one.
            if max_dist <= radius_bound {
                continue;
            }
            let ok = region.contains(cx, cy);
            if ok && d > best_dist {
                best_x = cx;
                best_y = cy;
                best_dist = d;
                found = true;
            }
            queue.push(GridCell {
                x: cx,
                y: cy,
                h_size: half,
                max_dist,
                support,
            });
        }
    }

    if !found {
        return None;
    }

    // Polish only when the search ran out of budget. A search that stopped on
    // the bound is already within `tolerance` of the optimum, so climbing it
    // again spends probes to recover nothing. A ridge outruns the budget before
    // the bisection closes and can sit a percent short.
    let (best_x, best_y) = if exhausted_budget {
        let (px, py, _) = polish(nb, region, best_x, best_y, best_dist, tolerance, h_size);
        (px, py)
    } else {
        (best_x, best_y)
    };

    let radius_sq = set.clearance_sq(best_x, best_y);
    if !radius_sq.is_finite() {
        return None;
    }

    Some(GridSolution {
        x: best_x,
        y: best_y,
        radius_sq: radius_sq.max(0.0),
        iterations,
        exhausted_budget,
    })
}

/// Climb the clearance field from an incumbent centre until no probe improves.
///
/// Each round probes the eight compass directions at a shrinking radius, and
/// accepts a probe only when it is inside the region and strictly better. The
/// radius halves whenever a round finds nothing, so the walk stops once the step
/// falls under `tol`.
fn polish(
    nb: &NearestObstacleIndex,
    region: &Region,
    start_x: f64,
    start_y: f64,
    start_d: f64,
    tol: f64,
    initial_step: f64,
) -> (f64, f64, f64) {
    let mut best_x = start_x;
    let mut best_y = start_y;
    let mut best_d = start_d;
    let mut step = initial_step;

    for _ in 0..256 {
        if step <= tol || !step.is_finite() {
            break;
        }
        let mut moved = false;
        for (dx, dy) in [
            (-1.0, 0.0),
            (1.0, 0.0),
            (0.0, -1.0),
            (0.0, 1.0),
            (-1.0, -1.0),
            (1.0, -1.0),
            (-1.0, 1.0),
            (1.0, 1.0),
        ] {
            let nx = best_x + dx * step;
            let ny = best_y + dy * step;
            let Some((d_sq, _)) = nb.nearest_distance_sq(nx, ny) else {
                continue;
            };
            let d = d_sq.max(0.0).sqrt();
            if d > best_d && region.contains(nx, ny) {
                best_x = nx;
                best_y = ny;
                best_d = d;
                moved = true;
            }
        }
        if !moved {
            step *= 0.5;
        }
    }
    (best_x, best_y, best_d)
}

#[inline]
fn nearest_distance(nb: &NearestObstacleIndex, x: f64, y: f64) -> f64 {
    match nb.nearest_distance_sq(x, y) {
        Some((d, _)) if d.is_finite() && d > 0.0 => d.sqrt(),
        _ => 0.0,
    }
}

/// Centroid of the obstacle set: the mean of every point and segment midpoint.
fn centroid_of(set: &ObstacleSet) -> Option<(f64, f64)> {
    let mut sx = 0.0;
    let mut sy = 0.0;
    let mut n = 0usize;
    for p in &set.points {
        sx += p[0];
        sy += p[1];
        n += 1;
    }
    for i in 0..set.segments.len() {
        sx += (set.segments.ax[i] + set.segments.bx[i]) * 0.5;
        sy += (set.segments.ay[i] + set.segments.by[i]) * 0.5;
        n += 1;
    }
    (n > 0).then(|| (sx / n as f64, sy / n as f64))
}
