<div align="center">

# MaxEmptyCircle

Finds the largest circle that fits between obstacles. The obstacles can be
points, lines, polygons, or any collection of them, and an optional boundary
restricts where the centre may sit.

[![crate](https://img.shields.io/crates/v/max-empty-circle.svg)](https://crates.io/crates/max-empty-circle)
[![docs](https://docs.rs/max-empty-circle/badge.svg)](https://docs.rs/max-empty-circle)
![MSRV](https://img.shields.io/badge/rustc-1.88+-ab6000.svg)
[![Tests](https://img.shields.io/badge/tests-55%20passing-brightgreen.svg)](https://github.com/Wolren/MaxEmptyCircle/actions/workflows/ci.yml)
[![GEOS](https://img.shields.io/badge/vs%20GEOS-5.04x%20faster-blue.svg)](https://github.com/Wolren/MaxEmptyCircle#performance)
[![License](https://img.shields.io/badge/license-Apache--2.0-blue.svg)](https://github.com/Wolren/MaxEmptyCircle?tab=License-1-ov-file)
[![Ko-fi](https://img.shields.io/badge/Ko--fi-Support%20Wolren-FF5E5B?logo=ko-fi&logoColor=white)](https://ko-fi.com/wolren)

</div>

A largest empty circle is the biggest circle that fits in open space. Every
obstacle limits it, so the circle touches the nearest one. This is the dual of
the maximum inscribed circle: that one sits inside a polygon and is limited by
the polygon's boundary, while this one sits between obstacles that need not
enclose the centre.

## API

| Function | Centre confined to |
|----------|--------------------|
| `largest_empty_circle` | Caller boundary, else the obstacle convex hull |
| `largest_empty_circle_slice` | Same, batch over rayon, input order preserved |

`MecResult` reports the centre, the radius and its square, the obstacles the
circle touches, the candidate count, and the solve path. `MecApproach` names the
path: `ExactVoronoi` or `GridSolver`.

Obstacles are any `geo_types::Geometry`: points, lines, polygons, multipolygons,
or a geometry collection whose parts are flattened into one set. The boundary
must be polygonal, because inside is undefined for a point or a line.

Obstacles that enclose no area, such as a single point or collinear points, give
a zero-radius circle centred on the obstacle. GEOS answers that case the same way.

Run `cargo run --example quickstart` for a complete program.

## How it solves

Point obstacles take an exact path. The maximum of a minimum of distances over a
closed region is attained at a Voronoi vertex, which is a Delaunay circumcentre,
or on the region boundary. One triangulation enumerates the interior candidates
in O(n log n); the boundary family adds ring corners, bisector-edge crossings
solved from the equidistance condition, and a per-edge refine within one sample
spacing. Each candidate is re-measured against every site, so a circumcentre with
a fourth site inside its circle scores lower and loses. The answer is exact, not
sampled.

Segments and polygons take a branch and bound over the region. The bound is the
1-Lipschitz property of obstacle distance: no point in a cell of half-size `h`
clears an obstacle by more than `d(c) + h*sqrt(2)`. Region membership is a filter,
not part of the bound, so a cell outside still holds the boundary optima. A
polygon contributes its edges rather than its interior, matching GEOS. When the
work cap trips, a rescan climbs from one representative per basin.

The reported radius is measured against the obstacle set at the returned centre,
never read off a bound.

## Correctness

`tests/geos_corpus.rs` ports the GEOS unit suite and its C API companion: 23 cases
carrying GEOS's own WKT and its published centres and radii. All 23 pass.

`tests/geos_parity.rs` cross-checks 10 cases against a live GEOS, asserting fit,
non-domination, and agreement.

`tests/corpus_invariants.rs` checks two properties. The radius must equal the
measured clearance at the reported centre, and no dense sample inside the region
may clear more than the returned circle. Its oracle shares no code with either
solver.

`tests/scaling.rs` fits the growth exponent over doubling input sizes and fails if
it reaches quadratic.

`tests/simd_equivalence.rs` pins the four-lane distance batch to the scalar result.

The workflow in `.github/workflows/ci.yml` runs the suite on every push and pull
request. It runs four feature sets, the 1.88 MSRV check, and a clippy, fmt, and
rustdoc job.

`examples/stress_parity.rs` runs the wide randomised sweep against GEOS:

| Sweep | Cases | Failures | Worst relative gap |
|-------|------:|---------:|-------------------:|
| Random point sets | 400 | 0 | 0 |
| Random segment sets | 200 | 0 | 1.3e-5 |
| Boundary-confined point sets | 100 | 0 | 0 |

A failure means the native radius was materially below GEOS's.

Three ported cases assert at least GEOS's radius rather than equality. GEOS's
published expectations for polygon obstacles are not optimal, because its cell
grid settles on a ridge before the true maximum. On 60 randomised polygon cases
the crate is better 17 times, equal 43, worse 0.

## Performance

The benchmark runs 14 cases against GEOS 3.14.1 through its C API. Only the solve
is timed, because the GEOS geometry is parsed before the clock starts. Both sides
get the same obstacles, region, and tolerance. Ratio is GEOS divided by native.
Both sides run serially, so the table measures the solver and not the threading.

| Case | Path | Primitives | native | GEOS | Ratio |
|------|------|-----------:|-------:|-----:|------:|
| points/8 | exact | 8 | 36 us | 37.0 ms | 1042x |
| points/32 | exact | 32 | 179 us | 60.6 ms | 338x |
| points/128 | exact | 128 | 2.53 ms | 17.6 ms | 7.0x |
| points/512 | grid | 512 | 17.4 ms | 30.5 ms | 1.8x |
| lattice/6x6 | exact | 36 | 138 us | 4.03 ms | 29x |
| lattice/12x12 | exact | 144 | 2.07 ms | 22.4 ms | 11x |
| lattice/20x20 | grid | 400 | 16.3 ms | 73.9 ms | 4.5x |
| buildings/3x3 | grid | 36 | 157 us | 870 us | 5.5x |
| buildings/6x6 | grid | 144 | 511 us | 4.88 ms | 9.6x |
| buildings/10x10 | grid | 400 | 1.55 ms | 20.1 ms | 13x |
| segments/16 | grid | 16 | 10.9 ms | 27.0 ms | 2.5x |
| segments/64 | grid | 64 | 11.7 ms | 70.6 ms | 6.0x |
| segments/256 | grid | 256 | 15.4 ms | 39.2 ms | 2.6x |
| confined/points64 | exact | 64 | 302 us | 305 us | 1.0x |

| Measure | Value |
|---------|------:|
| total, 14 cases | 5.04x |
| cases native solves faster | 14 of 14 |
| cases native radius materially worse | 0 |

The ratio varies by case. Across the 14:

| Per-case ratio | Value |
|----------------|------:|
| p10 | 1.79x |
| median | 7.73x |
| p90 | 309x |

The point-heavy rows are flattered, because GEOS pays a large fixed cost in its
branch-and-bound: its 37 ms for 8 sites is mostly overhead. The band on realistic
work is 1.8x to 13x.

Growth stays linear, fitted as `t = c * n^k`:

| Sweep | Range | Exponent |
|-------|-------|---------:|
| Point obstacles | 10 to 100k | k = 1.05 |
| Segment obstacles | 10 to 10k | k = 0.20 |

Segment time is dominated by the recovery polish, which is bounded by the region
rather than the obstacle count, so it is flat until the index cost takes over.

The two solvers agree on the answer. Across the 14 cases the largest radius gap is
8.6e-8 relative, and across the 700 randomised cases it is 1.3e-5.

Regenerate with
`cargo run --release --features bench-geos --example compare_geos`. The build
needs the GEOS C library. `build.rs` reads `GEOS_LIB_DIR`, or falls back to
`geos-config` on `PATH`.

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

The `python` feature builds the extension module `max_empty_circle`. Build it with
maturin, or with `cargo build --features python,extension-module`.

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
| `geo-traits` | Generic solves over any `geo_traits` geometry | no |
| `simd` | Four-lane `wide` batch for the climb probes | no |
| `shewchuk` | Exact `orient2d` predicates from `robust` in containment | no |
| `hotpath` | `hotpath::measure` instrumentation on the solver entry points | no |
| `python` | PyO3 extension module named `max_empty_circle` | no |
| `extension-module` | Build the Python module as an importable extension | no |
| `bench-geos` | Comparison only: links the GEOS C library for benchmarks and parity tests | no |

## License

Apache-2.0
