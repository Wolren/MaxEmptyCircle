<div align="center">

# MaxEmptyCircle

Finds the largest circle that fits between obstacles. The obstacles can be
points, lines, polygons, or any collection of them, and an optional boundary
restricts where the centre may sit.

[![MSRV](https://img.shields.io/badge/rustc-1.88+-ab6000.svg)]()
[![License](https://img.shields.io/badge/license-Apache--2.0-blue.svg)](LICENSE)

</div>

A largest empty circle is the biggest circle that fits in open space. Every
obstacle limits it, so the circle touches the nearest one. This is the dual of
the maximum inscribed circle: that one sits inside a polygon and is limited by
the polygon's own boundary, while this one sits between obstacles that need not
enclose the centre at all.

## API

| Function | Use |
|----------|-----|
| `largest_empty_circle` | one obstacle geometry, optional boundary |
| `largest_empty_circle_slice` | rayon-parallel batch, input order preserved |

`MecResult` reports the centre, the radius and its square, the obstacles the
circle touches, the candidate count, and the solve path. `MecApproach` names the
path: `ExactVoronoi` or `GridSolver`.

Obstacles are any `geo_types::Geometry`: points, lines, polygons, multipolygons,
or a geometry collection whose parts are flattened into one obstacle set. The
boundary must be polygonal, because "inside" is undefined for a point or a line.

Run `cargo run --example quickstart` for a complete program.

## Why a boundary is not optional in spirit

An empty circle between two finite obstacles is unbounded. Two parallel walls
have an arbitrarily large empty circle far from both, so an unconstrained centre
has no maximum, and the problem has no answer.

So the region matters, and the crate makes the rule explicit rather than silent.
With a caller boundary the centre is confined to it. Without one the centre is
confined to the convex hull of the obstacles, which is the GEOS default: the
answer then lies among the obstacles, which is what a facility-siting or
label-placement question usually means. Obstacles that enclose no area, such as a
single point or a set of collinear points, are rejected with an error naming the
reason, because no finite answer exists for them.

## How it solves

The crate picks a path from the obstacle types.

**Exact Voronoi, for point obstacles.** The largest empty circle among point
sites has its centre at a Voronoi vertex, so the answer is found by enumerating
circumcentres of site triples and re-measuring each candidate against every site.
The measurement is what makes it exact: a circumcentre with a fourth site inside
its circle scores lower and loses, so degeneracies filter themselves out. Pair
midpoints cover the two-site case. This path is exact, not approximated.

Confinement adds one candidate family that a vertex enumeration alone misses. The
clearance field increases outward across the region boundary, so a constrained
optimum can sit exactly on it, where no Voronoi vertex reaches. The crate solves
the site-pair bisector against each region edge analytically for those, and sweeps
the region as a safety net.

**Grid search, for segments and polygons.** A branch and bound over the region,
driven by the 1-Lipschitz property of obstacle distance: no point in a cell of
half-size `h` can clear an obstacle by more than `d(c) + h*sqrt(2)`. A compass
hill climb finishes the search when the cell budget runs out first.

The reported radius is always measured against the obstacle set at the returned
centre, never inferred from a bound, so a returned circle is empty by
construction.

## Correctness

`tests/corpus_invariants.rs` checks three properties, two of which found real
bugs during development.

Fit: the reported radius must equal the measured clearance at the reported
centre. Reach: no dense sample inside the region may clear more than the returned
circle. Route: point-only input must take the exact path, and the exact answer
must not lose to an independent brute-force sample.

The oracle in that file is separate from the solver: it samples the obstacle hull
on a dense grid and climbs, sharing no code with either path. It caught a missing
boundary-candidate family, and it caught a point-in-region test that rejected the
region boundary itself.

## Features

| Feature | Description | Default |
|---------|-------------|---------|
| `geo-traits` | generic solves over any `geo_traits` geometry | no |
| `shewchuk` | exact `orient2d` predicates from `robust` in the containment test | no |
| `hotpath` | `hotpath::measure` instrumentation on the solver entry points | no |
| `python` | PyO3 extension module named `max_empty_circle` | no |
| `extension-module` | build the Python module as an importable extension | no |
| `bench-geos` | comparison only: links the GEOS C library so benchmarks can cross-check | no |

## Python bindings

The `python` feature builds the extension module `max_empty_circle`. Build it
with maturin, or with `cargo build --features python,extension-module`.

| Function | Signature |
|----------|-----------|
| `solve_mec_py` | `(points, segments=None, boundary=None)` |

Points are `(x, y)` tuples, segments are `(x0, y0, x1, y1)` tuples, and the
boundary is a closed ring of `(x, y)` tuples. A rejected input raises
`ValueError`.

## C API

`src/ffi.rs` is always compiled. It exposes `mec_solve` and `mec_options_default`.
Flat `[x, y, ...]` point and `[x0, y0, x1, y1, ...]` segment arrays go in, and
`#[repr(C)]` structs come out. Each function returns 0 on success and -1 on
failure, a panic becomes -1, and a null or short buffer is rejected instead of
crashing.

## License

Apache-2.0
