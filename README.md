<div align="center">

# MaxEmptyCircle

Finds the largest circle that fits between obstacles. The obstacles can be
points, lines, polygons, or any collection of them, and an optional boundary
restricts where the centre may sit.

[![crate](https://img.shields.io/crates/v/max-empty-circle.svg)](https://crates.io/crates/max-empty-circle)
[![docs](https://docs.rs/max-empty-circle/badge.svg)](https://docs.rs/max-empty-circle)
![MSRV](https://img.shields.io/badge/rustc-1.88+-ab6000.svg)
[![Tests](https://img.shields.io/badge/tests-28%20passing-brightgreen.svg)](https://github.com/Wolren/MaxEmptyCircle/actions/workflows/ci.yml)
[![GEOS](https://img.shields.io/badge/vs%20GEOS-4.08x%20faster-blue.svg)](https://github.com/Wolren/MaxEmptyCircle#performance)
[![License](https://img.shields.io/badge/license-Apache--2.0-blue.svg)](https://github.com/Wolren/MaxEmptyCircle?tab=License-1-ov-file)
[![Ko-fi](https://img.shields.io/badge/Ko--fi-Support%20Wolren-FF5E5B?logo=ko-fi&logoColor=white)](https://ko-fi.com/wolren)

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
crossings computed from the equidistance condition, or the maximum of one distance
function along an edge, refined within one sample spacing. Every candidate is
re-measured against every site, so a circumcentre with a fourth site inside its
circle scores lower and loses, which is what handles lattice degeneracies. The
module documentation carries the full completeness argument. This path is exact,
not approximated.

**Grid search, for segments and polygons.** A branch and bound over the region,
driven by the 1-Lipschitz property of obstacle distance: no point in a cell of
half-size `h` can clear an obstacle by more than `d(c) + h*sqrt(2)`. Region
membership is a filter rather than part of the bound, so a cell outside still
carries a valid bound and still holds the boundary optima. When the work cap
trips, a coarse rescan of the region climbs from one representative of each basin,
which is what stops a single walk settling in the wrong one.

The reported radius is always measured against the obstacle set at the returned
centre, never inferred from a bound, so a returned circle is empty by
construction.

## Correctness

Two suites, both runnable locally.

`tests/corpus_invariants.rs` checks fit, meaning the reported radius equals the
measured clearance at the reported centre, and reach, meaning no dense sample
inside the region clears more than the returned circle. Its oracle is a coarse
grid plus a climb, sharing no code with either solver, and it caught two real
defects during development: a search that refined outside the region forever, and
a boundary optimum that a vertex enumeration cannot express.

`tests/geos_parity.rs` runs 10 cases against GEOS and asserts fit,
non-domination, and agreement on each. A wider sweep is what the numbers below
come from, and it is what found two defects the narrow suite missed.

`examples/stress_parity.rs` reproduces that sweep:

| Sweep | Cases | Failures | Worst relative gap |
|-------|------:|---------:|-------------------:|
| Random point sets | 400 | 0 | 0 |
| Random segment sets | 200 | 0 | 1.3e-5 |
| Boundary-confined point sets | 100 | 0 | 0 |

A failure means the native radius was materially below GEOS's. The worst segment
gap is 1.3e-5 relative, inside the tolerance the grid resolves to.

## Performance

Measured against GEOS 3.14.1 through its C API, on the same obstacles, the same
region and the same tolerance. Only the solve is timed. Every row reports the
radius both sides returned, so a speed claim is only made where the answers agree.

| Case | Path | Primitives | native | GEOS | Ratio |
|------|------|-----------:|-------:|-----:|------:|
| points/8 | exact | 8 | 32 us | 34.5 ms | 1071x |
| points/32 | exact | 32 | 174 us | 59.0 ms | 339x |
| points/128 | exact | 128 | 2.01 ms | 15.7 ms | 7.8x |
| points/512 | grid | 512 | 27.4 ms | 30.2 ms | 1.1x |
| lattice/6x6 | exact | 36 | 122 us | 3.88 ms | 32x |
| lattice/12x12 | exact | 144 | 2.02 ms | 20.9 ms | 10x |
| lattice/20x20 | grid | 400 | 4.53 ms | 68.9 ms | 15x |
| buildings/3x3 | grid | 36 | 327 us | 852 us | 2.6x |
| buildings/6x6 | grid | 144 | 548 us | 4.91 ms | 9.0x |
| buildings/10x10 | grid | 400 | 1.52 ms | 35.1 ms | 23x |
| segments/16 | grid | 16 | 17.5 ms | 25.3 ms | 1.5x |
| segments/64 | grid | 64 | 18.6 ms | 66.7 ms | 3.6x |
| segments/256 | grid | 256 | 24.0 ms | 38.3 ms | 1.6x |
| confined/points64 | exact | 64 | 295 us | 315 us | 1.1x |

**Total 4.08x, median 8.96x, 14 of 14 cases faster, no case materially worse
than GEOS.**

The point-heavy rows are flattered: GEOS pays a large fixed cost in its
branch-and-bound, so its 34 ms for 8 sites is mostly overhead. The honest reading
of the headline is the 1.5x to 23x band on realistic work.

Regenerate with:

```
GEOS_LIB_DIR=D:\Miniconda\Library\lib \
  cargo run --release --features bench-geos --example compare_geos
```

The build needs the GEOS C library. `build.rs` finds it through `GEOS_LIB_DIR`,
or through `geos-config` on `PATH`.

## Integration with the geo ecosystem

Enable `geo-traits` to accept a foreign geometry directly. WKB views, GeoArrow
arrays, and GeoZero streams then need no conversion.

| Function | Accepts |
|----------|---------|
| `largest_empty_circle_trait` | any `GeometryTrait` with `f64` coordinates |
| `largest_empty_circle_trait_with_boundary` | plus any `PolygonTrait` boundary |

A foreign geometry with no `geo_types` equivalent, such as an empty point, is
rejected with an error naming the reason rather than panicking.

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
`#[repr(C)]` structs come out.

Each function returns 0 on success and -1 on failure. A panic becomes -1, and a
null or short buffer is rejected instead of crashing. `MecOptionsC` carries no
tunable field yet. It exists so the ABI has a stable options pointer.

## Features

| Feature | Description | Default |
|---------|-------------|---------|
| `bench-geos` | Comparison only: links the GEOS C library so benchmarks and the parity suite can cross-check the native solver. Never compiled into the solve path. | no |
| `geo-traits` | Generic solves over any `geo_traits` geometry | no |
| `shewchuk` | Exact `orient2d` predicates from `robust` in the containment test | no |
| `hotpath` | `hotpath::measure` instrumentation on the solver entry points | no |
| `python` | PyO3 extension module named `max_empty_circle` | no |
| `extension-module` | Build the Python module as an importable extension | no |

## License

Apache-2.0
