//! The two solve paths, plus the GEOS comparison behind `bench-geos`.

pub mod apex;
#[cfg(feature = "bench-geos")]
pub mod bench_geos;
pub mod exact;
pub mod grid;
#[cfg(feature = "simd")]
pub mod simd_probe;
