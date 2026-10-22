//! The two solve paths, plus the GEOS comparison behind `bench-geos`.

#[cfg(feature = "bench-geos")]
pub mod bench_geos;
pub mod exact;
pub mod grid;
