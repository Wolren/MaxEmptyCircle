<div align="center">

# MaxEmptyCircle

Finds the largest circle that fits between obstacles.

![MSRV](https://img.shields.io/badge/rustc-1.88+-ab6000.svg)
[![Tests](https://img.shields.io/badge/tests-57%20passing-brightgreen.svg)](https://github.com/Wolren/MaxEmptyCircle/actions/workflows/ci.yml)
[![GEOS](https://img.shields.io/badge/vs%20GEOS-6.32x%20faster-blue.svg)](https://github.com/Wolren/MaxEmptyCircle#performance)
[![License](https://img.shields.io/badge/license-Apache--2.0-blue.svg)](https://github.com/Wolren/MaxEmptyCircle?tab=License-1-ov-file)
[![Ko-fi](https://img.shields.io/badge/Ko--fi-Support%20Wolren-FF5E5B?logo=ko-fi&logoColor=white)](https://ko-fi.com/wolren)

</div>

The largest empty circle is the biggest circle that fits in open space. Obstacles
are anything the circle must avoid: points, lines, or polygons, alone or mixed. An
optional boundary restricts where the centre may sit, and it is a separate input
from the obstacles: a building is an obstacle, a property line is a boundary. The
reported radius is the measured clearance at the reported centre, so a returned
circle is empty by construction. Input that encloses no area gives a zero-radius
circle on the obstacle, as GEOS does.

## Correctness

`tests/geos_corpus.rs` ports the GEOS unit suite and its C API companion: 23 cases
carrying GEOS's own WKT with its published centres and radii. All 23 pass.

Three of those assert at least GEOS's radius rather than equality. GEOS's published
expectations for polygon obstacles are not optimal, because its cell grid settles
on a ridge before the true maximum. On 60 randomised polygon cases this crate is
better 17 times, equal 43, worse 0.

`tests/corpus_invariants.rs` checks two properties. The radius must equal the
measured clearance at the reported centre, and no dense sample inside the region
may clear more than the returned circle. Its oracle shares no code with either
solver.

`tests/geos_parity.rs` cross-checks 10 cases against a live GEOS, asserting fit,
non-domination, and agreement.

A randomised sweep of 900 cases runs in `examples/stress_parity.rs`:

| Sweep | Cases | Failures | Worst relative gap |
|-------|------:|---------:|-------------------:|
| Random point sets | 400 | 0 | 0 |
| Random segment sets | 200 | 0 | 0 |
| Boundary-confined point sets | 100 | 0 | 0 |
| Mixed points and segments | 200 | 0 | 0 |

A failure means the native radius was materially below GEOS's, measured at 1e-6
relative, and the worst relative gap is the largest amount by which GEOS beat the
native answer. Both are zero on every sweep: no case loses to GEOS, and in the
cases where the two disagree this crate returns the larger circle. One-off runs of
10,000 segment and polygon cases and 4,000 sparse and mixed cases reported the same.

The gap used to be a per-mille on a handful of segment sets, and it was structural
rather than a tuning miss. A maximum of the clearance field is either equidistant
from three obstacles or pinned to the region boundary, and the grid path used to
approach both with a walk. A walk cannot resolve either one: on a peak the field is
flat to first order, and at a boundary optimum the only improving direction runs
exactly along the region edge, which a compass step cannot take. Both families are
now enumerated in closed form instead, which is why the sweep is clean rather than
merely inside a tolerance.

## Performance

i5-12400F (6C/12T), release, MSVC, against GEOS 3.13.1 through its C API. Only the
solve is timed: the GEOS geometry is parsed before the clock starts. Both sides get
the same obstacles, region, and tolerance. Ratio is GEOS divided by native, so above
1 means this crate wins. Both sides run serially, so the table measures the solver
and not the threading.

Every figure is the best of five runs of a harness that itself takes the best of
three repetitions, because a timing only drifts upward: the minimum is the cleanest
estimate of the solve. Single-run figures published earlier are correspondingly
higher, on both sides.

| Case | Path | Primitives | MaxEmptyCircle | GEOS | Ratio |
|------|------|-----------:|---------------:|-----:|------:|
| points/8 | exact | 8 | 31 us | 33.4 ms | 1087x |
| points/32 | exact | 32 | 167 us | 56.8 ms | 340x |
| points/128 | exact | 128 | 1.91 ms | 15.9 ms | 8.35x |
| points/512 | grid | 512 | 12.3 ms | 28.4 ms | 2.31x |
| lattice/6x6 | exact | 36 | 122 us | 3.90 ms | 32x |
| lattice/12x12 | exact | 144 | 1.91 ms | 21.3 ms | 11.2x |
| lattice/20x20 | grid | 400 | 14.0 ms | 68.4 ms | 4.87x |
| buildings/3x3 | grid | 36 | 531 us | 782 us | 1.47x |
| buildings/6x6 | grid | 144 | 619 us | 4.49 ms | 7.25x |
| buildings/10x10 | grid | 400 | 1.64 ms | 17.9 ms | 10.9x |
| segments/16 | grid | 16 | 7.49 ms | 25.0 ms | 3.34x |
| segments/64 | grid | 64 | 8.47 ms | 67.6 ms | 7.98x |
| segments/256 | grid | 256 | 10.7 ms | 36.9 ms | 3.44x |
| confined/points64 | exact | 64 | 280 us | 210 us | 0.75x |

Total 6.32x, 13 of 14 cases faster, no case materially worse. The ratio varies by
case: p10 1.47x, median 7.98x, p90 1087x. The spread is wider than it looks, because
the ratio is a speed-up and the slowest case is bounded below by 1.

The three `buildings/*` rows carry the exact refinement as a fixed cost per solve,
so they sit a few percent above an earlier table's figures for the same shapes. The
segment and large point-set rows sit below it, because the refinement replaced a
per-basin descent that cost more than the refinement does.

`confined/points64` is the one case where GEOS wins, at 0.75x. It is a 64-site solve
inside a small boundary, so the exact path runs its full boundary family against a
short ring and GEOS's grid finishes first. An earlier run reported it at 1.3x because
the timed GEOS region included building and parsing the WKT, which is setup rather
than solver work; the geometry is now prepared before the clock starts, so the figure
is the solve alone.

The point-heavy rows are flattered, because GEOS pays a large fixed cost in its
branch-and-bound: its 33 ms for 8 sites is mostly overhead. The band on realistic
work is 1.5x to 11x.

Growth stays linear, fitted as `t = c * n^k` over doubling input sizes:

| Sweep | Range | Exponent |
|-------|-------|---------:|
| Point obstacles | 10 to 100k | k = 1.03 |
| Segment obstacles | 10 to 10k | k = 0.17 |

`tests/scaling.rs` asserts both exponents stay under the subquadratic bound, so a
complexity regression fails CI. Segment time is dominated by the recovery rescan,
which is bounded by the region rather than the obstacle count, so it is flat until
the index cost takes over.

The two solvers agree on the answer. The largest radius gap is 7.2e-8 relative
across the 14 cases, and no randomised case out of the 900 is below GEOS at all. The
point path is exact, so agreement there shows the exact enumeration and the
iterative search reach the same optimum. The segment path is now enumerated for both
families of maxima, so agreement means GEOS's own tolerance rather than a shared
heuristic: where the two disagree, this crate returns the larger circle.

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
boundary is a closed ring of `(x, y)` tuples. The result reports the centre, the
radius, and the path. A rejected input raises `ValueError`.

## C API

`src/ffi.rs` is always compiled, so the C ABI needs no feature. Build with
`cargo build --release`.

| Symbol | Signature |
|--------|-----------|
| `mec_solve` | `(points, points_len, segments, segments_len, boundary, boundary_len, options, out) -> c_int` |
| `mec_options_default` | `() -> MecOptionsC` |

Flat `[x, y, ...]` point and `[x0, y0, x1, y1, ...]` segment arrays go in, and a
zero `boundary_len` means no boundary.
`MecResultC` and `MecOptionsC` are `#[repr(C)]`. Each function returns 0 on
success and -1 on failure. A panic becomes -1, and a null or short buffer is
rejected instead of crashing.

## Features

| Feature | Description | Default |
|---------|-------------|---------|
| `geo-traits` | Generic solves over any `geo_traits` geometry | no |
| `simd` | Four-lane `wide` batch for the climb probes | no |
| `shewchuk` | Exact `orient2d` predicates from `robust` in the containment test | no |
| `hotpath` | `hotpath::measure` instrumentation on the solver entry points | no |
| `python` | PyO3 extension module named `max_empty_circle` | no |
| `extension-module` | Build the Python module as an importable extension | no |
| `bench-geos` | Comparison only: links the GEOS C library for benchmarks and parity tests | no |

`rayon` and `libc` are unconditional dependencies. The parallel batch
`largest_empty_circle_slice` and the C ABI need no feature.

## License

Apache-2.0
