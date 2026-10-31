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

The largest empty circle is the biggest circle that fits in open space. Obstacles
are anything the circle must avoid (points, lines, polygons, alone or mixed), and an
optional boundary confines the centre. The reported radius is the measured clearance
at the returned centre, so a returned circle is empty by construction. Input that
encloses no area gives a zero-radius circle on the obstacle, as GEOS does. The API
reference and the algorithm detail live in the crate documentation, under
`cargo doc --open`.

## Correctness

`tests/corpus_invariants.rs` checks fit and reach. The reported radius must equal the
measured clearance at the reported centre, and no dense sample inside the region may
clear more than the returned circle. Its oracle shares no code with either solver.

`tests/geos_corpus.rs` ports the GEOS unit suite, 23 cases on GEOS's own WKT. All 23
pass, and three assert at least GEOS's radius rather than equality, because the
published expectations for polygon obstacles are not optimal.

`tests/geos_parity.rs` cross-checks 10 cases against a live GEOS.
`examples/stress_parity.rs` sweeps 900 randomised cases.

| Sweep | Cases | Failures | Worst gap |
|-------|------:|---------:|----------:|
| Points | 400 | 0 | 0 |
| Segments | 200 | 0 | 0 |
| Confined points | 100 | 0 | 0 |
| Points and segments | 200 | 0 | 0 |

Zero failures at 1e-6 relative, and a zero worst gap, so wherever the two answers
differ this crate returns the larger circle. One-off runs of 10,000 segment and
polygon cases and 4,000 mixed cases agree.

CI runs four feature sets plus the 1.88 MSRV check and clippy, fmt, and rustdoc. No
check needs an external solver.

## Performance

14 cases in `examples/compare_geos.rs`, covering both solve paths, with the same
obstacles, region, and tolerance for each solver. Only the solve is timed. Measured
on an i5-12400F, release, MSVC, against GEOS 3.13.1.

| Measure | MaxEmptyCircle | GEOS | Ratio |
|---------|---------------:|-----:|------:|
| total | **60.3 ms** | 381.0 ms | **6.32x** |
| cases native solves faster | 13 of 14 | 1 of 14 | |

Per-case ratio: p10 1.47x, median 7.98x, p90 1087x. The point rows dominate that
spread because GEOS pays a fixed cost per solve.

| Case | Path | Prims | MaxEmptyCircle | GEOS | Ratio |
|------|------|------:|---------------:|-----:|------:|
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

The exact refinement costs a fixed amount per solve, so `buildings/3x3` is the
closest grid row at 1.47x. `confined/points64` is the one case GEOS wins, at 0.75x,
because a 64-site solve inside a small boundary sends the exact path through its
whole boundary family against a short ring.

Figures are the best of five runs, each already the best of three repetitions.
Repeat runs land between 6.32x and 6.33x.

Growth is linear, `k = 1.03` for points over 10 to 100k and `k = 0.17` for segments
over 10 to 10k, both under the bound `tests/scaling.rs` enforces. Across the 14 cases
the largest radius gap is 7.2e-8 relative, and no randomised case lands below GEOS.

Regenerate with `cargo run --release --features bench-geos --example compare_geos`,
against a GEOS C library found through `GEOS_LIB_DIR` or `geos-config`.

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

Points are `(x, y)` tuples, segments are `(x0, y0, x1, y1)` tuples, and the boundary
is a closed ring of `(x, y)` tuples. A rejected input raises `ValueError`.

## C API

`src/ffi.rs` is always compiled, so the C ABI needs no feature.

| Symbol | Signature |
|--------|-----------|
| `mec_solve` | `(points, points_len, segments, segments_len, boundary, boundary_len, options, out) -> c_int` |
| `mec_options_default` | `() -> MecOptionsC` |

Flat `[x, y, ...]` point and `[x0, y0, x1, y1, ...]` segment arrays go in, and a zero
`boundary_len` means no boundary. `MecResultC` and `MecOptionsC` are `#[repr(C)]`.
Each function returns 0 on success and -1 on failure, and a panic becomes -1.

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

`rayon` and `libc` are unconditional dependencies.

## License

Apache-2.0
