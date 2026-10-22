//! Build script for the `bench-geos` feature.
//!
//! Only active when the feature is on. Links the GEOS C library so the
//! comparison harness can call `GEOSLargestEmptyCircle_r`.
//!
//! Discovery order: `GEOS_LIB_DIR`, then `geos-config` on `PATH`. `pkg-config`
//! is not used, because `geos-sys` is not a dependency and this crate needs only
//! the link line, not version detection.

fn main() {
    println!("cargo:rerun-if-changed=build.rs");
    println!("cargo:rerun-if-env-changed=GEOS_LIB_DIR");

    if std::env::var_os("CARGO_FEATURE_BENCH_GEOS").is_none() {
        return;
    }

    if let Some(dir) = std::env::var_os("GEOS_LIB_DIR") {
        let dir = dir.to_string_lossy().to_string();
        println!("cargo:rustc-link-search=native={dir}");
        println!("cargo:rustc-link-lib=dylib=geos_c");
        return;
    }

    // Fall back to geos-config, which most installs ship.
    match std::process::Command::new("geos-config")
        .arg("--libs")
        .output()
    {
        Ok(out) if out.status.success() => {
            let libs = String::from_utf8_lossy(&out.stdout);
            for token in libs.split_whitespace() {
                if let Some(path) = token.strip_prefix("-L") {
                    println!("cargo:rustc-link-search=native={path}");
                }
            }
            println!("cargo:rustc-link-lib=dylib=geos_c");
        }
        _ => {
            println!(
                "cargo:warning=bench-geos is on but GEOS was not found. \
                 Set GEOS_LIB_DIR to the directory holding geos_c, or put \
                 geos-config on PATH. The feature will fail to link."
            );
            println!("cargo:rustc-link-lib=dylib=geos_c");
        }
    }
}
