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

**Exact Voronoi, for point obstacles.** The maximum of a minimum of distances over
a closed region is attained either at an interior point, where three sites must be
equidistant and none closer, or on the region boundary. Interior optima are
Voronoi vertices, which are the Delaunay circumcentres, so one triangulation at
`O(n log n)` enumerates them all. Boundary optima are corners, bisector-edge
crossings, or the maximum of one distance function along an edge, and each is
computed directly rather than sampled. Every candidate is re-measured against
every site, so a circumcentre with a fourth site inside its circle scores lower
and loses, which is what handles lattice degeneracies. The module documentation
carries the full completeness argument. This path is exact, not approximated.

**Grid search, for segments and polygons.** A branch and bound over the region,
driven by the 1-Lipschitz property of obstacle distance: no point in a cell of
half-size `h` can clear an obstacle by more than `d(c) + h*sqrt(2)`. A compass
hill climb finishes the search when the cell budget runs out first.

The reported radius is always measured against the obstacle set at the returned
centre, never inferred from a bound, so a returned circle is empty by
construction.

## Performance

Measured against GEOS 3.14.1 through its C API, on the same obstacles, the same
region and the same tolerance. Only the solve is timed. Every row reports the
radius both sides returned, so a speed win is only claimed where the answers
agree.

| Case | Path | Primitives | native | GEOS | Ratio |
|------|------|-----------:|-------:|-----:|------:|
| points/8 | exact | 8 | 32 us | 34.3 ms | 1064x |
| points/32 | exact | 32 | 172 us | 56.1 ms | 326x |
| points/128 | exact | 128 | 1.89 ms | 16.5 ms | 8.7x |
| points/512 | grid | 512 | 14.5 ms | 29.9 ms | 2.1x |
| lattice/6x6 | exact | 36 | 108 us | 3.94 ms | 36x |
| lattice/12x12 | exact | 144 | 1.0 ms | 21.0 ms | 21x |
| lattice/20x20 | grid | 400 | 4.54 ms | 68.2 ms | 15x |
| buildings/3x3 | grid | 36 | 316 us | 859 us | 2.7x |
| buildings/6x6 | grid | 144 | 541 us | 4.84 ms | 9.0x |
| buildings/10x10 | grid | 400 | 1.55 ms | 19.2 ms | 12x |
| segments/16 | grid | 16 | 9.88 ms | 25.7 ms | 2.6x |
| segments/64 | grid | 64 | 10.0 ms | 69.2 ms | 6.9x |
| segments/256 | grid | 256 | 14.8 ms | 36.8 ms | 2.5x |
| confined/points64 | exact | 64 | 267 us | 306 us | 1.1x |

**Total 6.5x, median 9.2x, 14 of 14 cases faster, zero radius discrepancies.**
Largest gap between the two radii across the suite is 8.6e-8 relative, which is
the tolerance regime rather than agreement by luck.

Regenerate with:

```
GEOS_LIB_DIR=D:\Miniconda\Library\lib \
  cargo run --release --features bench-geos --example compare_geos
```

The build needs the GEOS C library. `build.rs` finds it through `GEOS_LIB_DIR`,
or through `geos-config` on `PATH`.

## Correctness

Two independent suites.

`tests/corpus_invariants.rs` checks that the reported radius equals the measured
clearance at the reported centre, and that no dense sample inside the region
clears more than the returned circle. Its oracle is a coarse grid plus a climb,
sharing no code with either solver, and it caught two real defects during
development: a search that refined outside the region forever, and a boundary
optimum that a vertex enumeration cannot express.

`tests/geos_parity.rs` runs 10 cases against GEOS, including 60 random point sets
and 30 random segment sets, and asserts three properties on each: fit, meaning the
radius is the measured clearance; non-domination, meaning the native radius is
never materially below GEOS's; and agreement within tolerance.

The exact path is not a heuristic, and the module documentation carries the
argument: the maximum of a minimum of distances over a closed region is attained
either at a Voronoi vertex, which the Delaunay triangulation enumerates, or on the
region boundary, where the field restricted to an edge is piecewise convex with
breakpoints at bisector crossings, so its maximum is at a corner, a crossing, or a
point found by refining within one sample spacing of the best sample.

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
