<div align="center">

# MaxEmptyCircle

Finds the largest circle that fits between obstacles.

[![crate](https://img.shields.io/crates/v/max-empty-circle.svg)](https://crates.io/crates/max-empty-circle)
[![docs](https://docs.rs/max-empty-circle/badge.svg)](https://docs.rs/max-empty-circle)
![MSRV](https://img.shields.io/badge/rustc-1.88+-ab6000.svg)
[![Tests](https://img.shields.io/badge/tests-57%20passing-brightgreen.svg)](https://github.com/Wolren/MaxEmptyCircle/actions/workflows/ci.yml)
[![GEOS](https://img.shields.io/badge/vs%20GEOS-6.32x%20faster-blue.svg)](https://github.com/Wolren/MaxEmptyCircle#performance)
[![License](https://img.shields.io/badge/license-Apache--2.0-blue.svg)](https://github.com/Wolren/MaxEmptyCircle?tab=License-1-ov-file)
[![Ko-fi](https://img.shields.io/badge/Ko--fi-Support%20Wolren-FF5E5B?logo=ko-fi&logoColor=white)](https://ko-fi.com/wolren)

</div>

The largest empty circle is the biggest circle that fits in open space, and
obstacles are anything it must avoid (points, lines, polygons, alone or mixed). An
optional boundary restricts where the centre may sit, and it is a separate input
from the obstacles, because a building is an obstacle while a property line is a
boundary. Since the reported radius is the measured clearance at the returned
centre, a returned circle is empty by construction. Input that encloses no area
yields a zero-radius circle on the obstacle, matching GEOS. The crate documentation
ships with the code and carries the API reference and the algorithm detail. Run
`cargo doc --open` and no feature is needed.

## Correctness

`tests/corpus_invariants.rs` checks two properties, and both come from the problem
rather than from the solver.

Fit comes first, which means the reported radius has to match the measured clearance
at the reported centre, even when that radius is the smaller of the two. Reach
comes second, and it forbids any dense sample inside the region from clearing more
than the returned circle. That oracle shares no code with either solver.

`tests/geos_corpus.rs` ports the GEOS unit suite and its C API companion, 23 cases
in all, each carrying GEOS's own WKT with its published centre and radius. All 23
pass. Three of them assert at least GEOS's radius rather than equality,
because the published expectations for polygon obstacles leave its grid short of
the true maximum.

`tests/geos_parity.rs` cross-checks 10 cases against a live GEOS, asserting fit,
non-domination, and agreement.

`examples/stress_parity.rs` sweeps 900 randomised cases.

| Sweep | Cases | Failures | Worst relative gap |
|-------|------:|---------:|-------------------:|
| Random point sets | 400 | 0 | 0 |
| Random segment sets | 200 | 0 | 0 |
| Boundary-confined point sets | 100 | 0 | 0 |
| Mixed points and segments | 200 | 0 | 0 |

A failure here means the native radius sat materially below GEOS's, at 1e-6
relative. Every sweep reports zero, so wherever the two answers part company this
crate returns the larger circle. One-off runs of 10,000 segment and polygon cases
and 4,000 sparse and mixed cases agree.

The grid path does not approach the maximum with a walk. A maximum of the clearance
field is either equidistant from three obstacles or pinned to the region boundary,
and a walk resolves neither: on a peak the field is flat to first order, and at a
boundary optimum the only improving direction runs along the region edge. Both
families are enumerated in closed form.

The workflow in `.github/workflows/ci.yml` runs the suite on every push and pull
request, covering four feature sets plus the 1.88 MSRV check and a clippy, fmt, and
rustdoc job. No CI check needs an external solver, because the GEOS comparison
stays behind its own feature.

## Performance

`examples/compare_geos.rs` holds 14 cases that cover both solve paths, and each
solver receives the same obstacles, region, and tolerance. Only the solve is timed,
since both sides build their geometry before the clock starts. Ratio is GEOS
divided by native, so above 1 means this crate wins, and both sides run serially to
measure the solver rather than the threading.

Figures below come from an i5-12400F (6C/12T), release, MSVC, against GEOS 3.13.1
through its C API.

| Measure | MaxEmptyCircle | GEOS | Ratio |
|---------|---------------:|-----:|------:|
| total, 14 cases | **60.3 ms** | 381.0 ms | **6.32x** |
| cases native solves faster | 13 of 14 | 1 of 14 | 93% |

| Per-case ratio | Value |
|----------------|------:|
| p10 | 1.47x |
| median | 7.98x |
| p90 | 1087x |

Each figure is the best of five runs of a runner that itself takes the best of
three repetitions, because timings only drift upward and the minimum is the
closest estimate of the solve. Repeat runs of this tree have landed between 6.32x
and 6.33x, and that spread makes the second digit noise.

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

The refinement costs a fixed amount per solve, so its share grows as the input
shrinks. buildings/3x3 is the smallest grid case and the closest row, at 1.47x.

`confined/points64` is the only case where GEOS wins at 0.75x. A 64-site solve
inside a small boundary sends the exact path through its full boundary family
against a short ring, while GEOS's grid finishes first. The timer covers the solve
alone, so the figure excludes the cost of building either geometry.

Point-heavy rows are flattered by a large fixed cost in GEOS's branch-and-bound,
and its 33 ms for 8 sites is mostly that overhead. Realistic work lands between
1.5x and 11x.

Growth stays linear, fitted as `t = c * n^k` over doubling input sizes.

| Sweep | Range | Exponent |
|-------|-------|---------:|
| Point obstacles | 10 to 100k | k = 1.03 |
| Segment obstacles | 10 to 10k | k = 0.17 |

`tests/scaling.rs` asserts both exponents stay under the subquadratic bound, so a
complexity regression fails CI. Segment time is dominated by the recovery rescan,
which is bounded by the region rather than the obstacle count, so it stays flat
until the index cost takes over.

The two solvers agree on the answer. The largest radius gap is 7.2e-8 relative
across the 14 cases, and no randomised case out of the 900 lands below GEOS.
On the point path that agreement shows the exact enumeration and the iterative
search reaching the same optimum. On the segment path it means GEOS's own tolerance
rather than a shared heuristic, because both families of maxima are enumerated.

Regenerate with
`cargo run --release --features bench-geos --example compare_geos`. That build needs
the GEOS C library, so `build.rs` reads `GEOS_LIB_DIR` or falls back to
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

`src/ffi.rs` is always compiled, so the C ABI needs no feature.

| Symbol | Signature |
|--------|-----------|
| `mec_solve` | `(points, points_len, segments, segments_len, boundary, boundary_len, options, out) -> c_int` |
| `mec_options_default` | `() -> MecOptionsC` |

Build with `cargo build --release`. Flat `[x, y, ...]` point and
`[x0, y0, x1, y1, ...]` segment arrays go in, and a zero `boundary_len` means no
boundary. `MecResultC` and `MecOptionsC` are `#[repr(C)]`. Each function returns 0
on success and -1 on failure. A panic becomes -1, and a null or short buffer is
rejected instead of crashing.

## Features

| Feature | Description | Default |
|---------|-------------|---------|
| `geo-traits` | Generic solves over any `geo_traits` geometry | no |
| `simd` | Four-lane `wide` batch for the climb probes | no |
| `shewchuk` | Exact `orient2d` predicates from the `robust` crate in the containment test | no |
| `hotpath` | `hotpath::measure` instrumentation on the solver entry points | no |
| `python` | PyO3 extension module named `max_empty_circle` | no |
| `extension-module` | Build the Python module as an importable extension | no |
| `bench-geos` | Comparison only. Links the GEOS C library for benchmarks and parity tests | no |

`rayon` and `libc` are unconditional dependencies. The parallel batch
`largest_empty_circle_slice` and the C ABI need no feature.

## License

Apache-2.0
