//! MEC tuning constants.
//!
//! Only the `MEC_*` family, mirroring the `MIC_*` layout of the sibling crate.

/// Half of the first-level cell count along the longer axis of the search box.
pub const MEC_GRID_SIDE: usize = 25;
/// Rounds the compass climb may take.
/// The walk also stops once it stops improving, so this only bounds a run whose
/// step keeps growing. A tight bound keeps a flat basin from paying for the whole
/// budget: the cost of one round is 32 distance queries.
pub const MEC_CLIMB_ROUNDS: usize = 1024;
/// Rounds of the shrink-only descent that refines the recovered peak.
pub const MEC_DESCENT_ROUNDS: usize = 512;
/// Consecutive rounds without improvement before the climb gives up.
pub const MEC_CLIMB_STAGNANT_ROUNDS: u32 = 64;
/// Grid cells popped before the work cap trips.
/// The cap only bounds a search that cannot converge on its bound. The recovery
/// polish supplies the accuracy, so a low cap costs nothing in correctness and
/// saves the whole population cost on the cases that trip it.
pub const MEC_GRID_CELL_BUDGET: usize = 5_000;
/// Relative gap below which a grid answer counts as equal to the exact seed.
/// The seed path reports candidates measured to this accuracy, so a grid that
/// lands inside the band has confirmed the candidate rather than beaten it.
pub const MEC_VORONOI_CONFIRM_REL_TOL: f64 = 1e-9;
/// Site count above which the exact triple enumeration declines.
/// The exact path is O(n^3) circumcentre solves each re-measured in O(n), so
/// this caps the work. Above it the grid path answers, which is approximate
/// but scales. 256 sites is ~2.8M triples, near the practical ceiling.
pub const MEC_EXACT_MAX_SITES: usize = 200;
/// Samples per axis of the region sweep the exact path runs before its climb.
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
/// The rescan only supplies starting basins; the climb supplies the accuracy, so
/// this stays coarse enough to cost less than the search it is recovering.
pub const MEC_POLISH_SCAN_SIDE: usize = 64;
/// Starting samples the polish climbs from.
pub const MEC_POLISH_SEEDS: usize = 40;
/// Divisions of the region span used as the basin separation radius.
/// Sampled points closer than `span / this` climb to the same peak, so the dedup
/// keeps one representative per neighbourhood instead of paying for every sample
/// inside it.
pub const MEC_POLISH_BASIN_DIVISIONS: usize = 48;
/// Relative band around the incumbent clearance whose primitives enter the exact
/// apex enumeration.
///
/// The apex is within one grid tolerance *in value* of the incumbent, but the two
/// can be far apart *in position* when the field is flat, and a primitive active at
/// the apex then sits further from the incumbent than the value gap suggests. The
/// band is set from the flat-field case: a peak that clears the incumbent by a
/// relative 1e-4 can have its third support up to about `sqrt(2 R delta) / R` away,
/// which for a shallow basin is a few per-mille. Two percent covers that with room,
/// and it costs nothing to be generous because the cap below, not the band, bounds
/// the work.
pub const MEC_APEX_BAND_REL: f64 = 2e-2;
/// Most primitives the apex enumeration takes, nearest first.
///
/// Each primitive contributes up to three carriers, and a triple can only use one
/// carrier per primitive, so the count is `27 * C(n, 3)`: six primitives is 540
/// triples, a few tens of microseconds. The cap is a ceiling on a degenerate tie,
/// not a normal case: a genuine apex is three-way, so a handful of primitives
/// covers it with spares.
pub const MEC_APEX_MAX_PRIMS: usize = 6;
