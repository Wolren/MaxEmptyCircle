//! MEC tuning constants.
//!
//! Only the `MEC_*` family, mirroring the `MIC_*` layout of the sibling crate.

/// Half of the first-level cell count along the longer axis of the search box.
pub const MEC_GRID_SIDE: usize = 25;
/// Cells the grid stage may pop before it returns its incumbent.
///
/// A ridge of equal clearance costs `length / tolerance` cells, so the budget
/// stops shapes the search cannot resolve and the hill climb finishes them.
pub const MEC_GRID_CELL_BUDGET: usize = 50_000;
/// Relative gap below which a grid answer counts as equal to the exact seed.
///
/// The seed path reports candidates measured to this accuracy, so a grid that
/// lands inside the band has confirmed the candidate rather than beaten it.
pub const MEC_VORONOI_CONFIRM_REL_TOL: f64 = 1e-9;
/// Site count above which the exact triple enumeration declines.
///
/// The exact path is O(n^3) circumcentre solves each re-measured in O(n), so
/// this caps the work. Above it the grid path answers, which is approximate
/// but scales. 256 sites is ~2.8M triples, near the practical ceiling.
pub const MEC_EXACT_MAX_SITES: usize = 200;
/// Samples per axis of the region sweep the exact path runs before its climb.
///
/// The constrained optimum can sit on the region boundary, where no Voronoi
/// vertex reaches, so the exact path sweeps the region and then climbs. This is
/// the coarse pass; the climb is what supplies the accuracy.
pub const MEC_EXACT_BOUNDARY_STEPS: usize = 192;
/// Hill-climb rounds the exact path spends after its region sweep.
pub const MEC_EXACT_CLIMB_ROUNDS: usize = 512;
/// Determinant magnitude below which a circumcentre solve counts as singular.
pub const MEC_EPS: f64 = 1e-12;
/// Samples per region edge before the golden-section refine.
pub const MEC_EDGE_SAMPLES: usize = 32;
/// Golden-section iterations refining each edge maximum.
pub const MEC_EDGE_REFINE: usize = 120;
/// Lattice samples per axis of the coarse rescan the polish runs.
///
/// The rescan only supplies starting basins; the climb supplies the accuracy, so
/// this stays coarse enough to cost less than the search it is recovering.
pub const MEC_POLISH_SCAN_SIDE: usize = 192;
/// Starting samples the polish climbs from.
pub const MEC_POLISH_SEEDS: usize = 40;
/// Divisions of the region span used as the basin separation radius.
///
/// Sampled points closer than `span / this` climb to the same peak, so the dedup
/// keeps one representative per neighbourhood instead of paying for every sample
/// inside it.
pub const MEC_POLISH_BASIN_DIVISIONS: usize = 48;
