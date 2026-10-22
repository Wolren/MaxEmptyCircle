//! Largest Empty Circle benchmark against GEOS, apples to apples.
//!
//! Both sides get the same obstacles, the same region, and the same tolerance,
//! and only the solve is timed. The native side runs the public entry point, so
//! its obstacle flatten is included; GEOS gets its WKT parsed before the clock
//! starts, because parsing is not what either solver is being measured on.
//!
//! Every row also reports the ratio of the two radii. A native win only means
//! something when the radii agree, so the summary counts any row where the
//! native radius is materially worse than GEOS's.
//!
//! Run:
//! ```text
//! GEOS_LIB_DIR=D:\Miniconda\Library\lib \
//!   cargo run --release --features bench-geos --example compare_geos
//! ```

use std::time::Instant;

use max_empty_circle::geo_types::{
    Coord, Geometry, LineString, MultiLineString, MultiPoint, Point, Polygon,
};
use max_empty_circle::mec::input::ObstacleSet;
use max_empty_circle::mec::solver::bench_geos::{geos_tolerance_for, solve_geos_with_tolerance};
use max_empty_circle::{largest_empty_circle, MecApproach, MecOptions};

/// A benchmark case: a name, its obstacles, and an optional boundary.
struct Case {
    name: String,
    obstacles: Geometry<f64>,
    boundary: Option<Geometry<f64>>,
}

struct Row {
    name: String,
    prims: usize,
    native_us: f64,
    geos_us: f64,
    native_r: f64,
    geos_r: f64,
    approach: MecApproach,
}

/// Deterministic xorshift so a run reproduces.
struct Rng(u64);

impl Rng {
    fn next(&mut self) -> u64 {
        self.0 ^= self.0 << 13;
        self.0 ^= self.0 >> 7;
        self.0 ^= self.0 << 17;
        self.0
    }

    fn unit(&mut self) -> f64 {
        (self.next() >> 11) as f64 / (1u64 << 53) as f64
    }

    fn range(&mut self, lo: f64, hi: f64) -> f64 {
        lo + (hi - lo) * self.unit()
    }
}

fn square(cx: f64, cy: f64, half: f64) -> Polygon<f64> {
    Polygon::new(
        LineString::from(vec![
            (cx - half, cy - half),
            (cx + half, cy - half),
            (cx + half, cy + half),
            (cx - half, cy + half),
            (cx - half, cy - half),
        ]),
        vec![],
    )
}

fn build_cases() -> Vec<Case> {
    let mut cases = Vec::new();

    // Point obstacles at several scales: the exact Voronoi path.
    for &n in &[8usize, 32, 128, 512] {
        let mut rng = Rng(0x9E37_79B9_7F4A_7C15 ^ n as u64);
        let pts: Vec<Point<f64>> = (0..n)
            .map(|_| Point::new(rng.range(0.0, 1000.0), rng.range(0.0, 1000.0)))
            .collect();
        cases.push(Case {
            name: format!("points/{n}"),
            obstacles: Geometry::MultiPoint(MultiPoint::new(pts)),
            boundary: None,
        });
    }

    // A lattice: every Voronoi vertex is degenerate, four sites equidistant.
    for &side in &[6usize, 12, 20] {
        let mut pts = Vec::new();
        for i in 0..side {
            for j in 0..side {
                pts.push(Point::new(i as f64 * 10.0, j as f64 * 10.0));
            }
        }
        cases.push(Case {
            name: format!("lattice/{side}x{side}"),
            obstacles: Geometry::MultiPoint(MultiPoint::new(pts)),
            boundary: None,
        });
    }

    // Building footprints: the realistic facility-siting case.
    for &grid_n in &[3usize, 6, 10] {
        let mut polys = Vec::new();
        for i in 0..grid_n {
            for j in 0..grid_n {
                polys.push(square(i as f64 * 60.0 + 10.0, j as f64 * 60.0 + 10.0, 14.0));
            }
        }
        cases.push(Case {
            name: format!("buildings/{grid_n}x{grid_n}"),
            obstacles: Geometry::MultiPolygon(geo_types::MultiPolygon::new(polys)),
            boundary: None,
        });
    }

    // Random segments: the case with no closed form at all.
    for &n in &[16usize, 64, 256] {
        let mut rng = Rng(0x2545_F491_4F6C_DD1D ^ n as u64);
        let segs: Vec<LineString<f64>> = (0..n)
            .map(|_| {
                let x = rng.range(0.0, 1000.0);
                let y = rng.range(0.0, 1000.0);
                LineString::from(vec![
                    Coord { x, y },
                    Coord {
                        x: x + rng.range(-80.0, 80.0),
                        y: y + rng.range(-80.0, 80.0),
                    },
                ])
            })
            .collect();
        cases.push(Case {
            name: format!("segments/{n}"),
            obstacles: Geometry::MultiLineString(MultiLineString::new(segs)),
            boundary: None,
        });
    }

    // Boundary-confined: the site must stay inside a parcel while clearing
    // obstacles that extend past it.
    {
        let mut rng = Rng(0xDEAD_BEEF_1234_5678);
        let pts: Vec<Point<f64>> = (0..64)
            .map(|_| Point::new(rng.range(0.0, 500.0), rng.range(0.0, 500.0)))
            .collect();
        cases.push(Case {
            name: "confined/points64".to_string(),
            obstacles: Geometry::MultiPoint(MultiPoint::new(pts)),
            boundary: Some(Geometry::Polygon(Polygon::new(
                LineString::from(vec![
                    (100.0, 100.0),
                    (260.0, 100.0),
                    (260.0, 260.0),
                    (100.0, 260.0),
                    (100.0, 100.0),
                ]),
                vec![],
            ))),
        });
    }

    cases
}

fn time_native(geom: &Geometry<f64>, boundary: Option<&Geometry<f64>>, reps: usize) -> (f64, f64) {
    let opts = MecOptions::default();
    let mut best = f64::INFINITY;
    let mut radius = 0.0;
    for _ in 0..reps {
        let start = Instant::now();
        let r = largest_empty_circle(geom, boundary, &opts).expect("native solve");
        let us = start.elapsed().as_secs_f64() * 1e6;
        if us < best {
            best = us;
        }
        radius = r.radius;
    }
    (best, radius)
}

fn time_geos(set: &ObstacleSet, boundary: Option<&Geometry<f64>>, reps: usize) -> (f64, f64) {
    let tolerance = geos_tolerance_for(set);
    let mut best = f64::INFINITY;
    let mut radius = 0.0;
    for _ in 0..reps {
        let start = Instant::now();
        let c = solve_geos_with_tolerance(set, boundary, tolerance).expect("geos solve");
        let us = start.elapsed().as_secs_f64() * 1e6;
        if us < best {
            best = us;
        }
        radius = c.radius;
    }
    (best, radius)
}

fn main() {
    let reps = 3;
    let cases = build_cases();
    let mut rows = Vec::new();

    for case in &cases {
        let set = match ObstacleSet::from_geometry(&case.obstacles) {
            Ok(s) => s,
            Err(e) => {
                eprintln!("skip {}: {e}", case.name);
                continue;
            }
        };
        let prims = set.point_count() + set.segments.len();
        let approach = match largest_empty_circle(
            &case.obstacles,
            case.boundary.as_ref(),
            &MecOptions::default(),
        ) {
            Ok(r) => r.approach,
            Err(e) => {
                eprintln!("skip {}: {e}", case.name);
                continue;
            }
        };
        let (native_us, native_r) = time_native(&case.obstacles, case.boundary.as_ref(), reps);
        let (geos_us, geos_r) = time_geos(&set, case.boundary.as_ref(), reps);
        rows.push(Row {
            name: case.name.clone(),
            prims,
            native_us,
            geos_us,
            native_r,
            geos_r,
            approach,
        });
    }

    println!();
    println!(
        "{:<22}{:>7}{:>9}{:>11}{:>11}{:>8}{:>10}{:>10}{:>9}",
        "case", "path", "prims", "native_us", "geos_us", "ratio", "native_r", "geos_r", "rel"
    );
    println!("{}", "-".repeat(97));

    let mut native_total = 0.0;
    let mut geos_total = 0.0;
    let mut worse = 0usize;
    let mut faster = 0usize;

    for row in &rows {
        let ratio = row.geos_us / row.native_us;
        let rel = (row.native_r - row.geos_r).abs() / row.geos_r.max(f64::MIN_POSITIVE);
        if row.native_us < row.geos_us {
            faster += 1;
        }
        // Materially worse: a gap far above either solver's tolerance.
        if row.native_r < row.geos_r * (1.0 - 1e-5) {
            worse += 1;
        }
        native_total += row.native_us;
        geos_total += row.geos_us;
        let path = match row.approach {
            MecApproach::ExactVoronoi => "exact",
            MecApproach::GridSolver => "grid",
            MecApproach::Geos => "geos",
        };
        println!(
            "{:<22}{:>7}{:>9}{:>11.1}{:>11.1}{:>7.2}x{:>10.3}{:>10.3}{:>8.1e}",
            row.name,
            path,
            row.prims,
            row.native_us,
            row.geos_us,
            ratio,
            row.native_r,
            row.geos_r,
            rel
        );
    }

    println!("{}", "-".repeat(97));
    println!(
        "{:<22}{:>7}{:>9}{:>11.1}{:>11.1}{:>7.2}x",
        "TOTAL",
        "",
        "",
        native_total,
        geos_total,
        geos_total / native_total
    );

    println!();
    println!("cases:                            {}", rows.len());
    println!("native faster:                    {faster}/{}", rows.len());
    println!("radius materially worse than GEOS: {worse}");
    println!(
        "exact-path cases:                 {}",
        rows.iter()
            .filter(|r| r.approach == MecApproach::ExactVoronoi)
            .count()
    );
    println!();

    let mut ratios: Vec<f64> = rows.iter().map(|r| r.geos_us / r.native_us).collect();
    ratios.sort_by(|a, b| a.partial_cmp(b).unwrap());
    if !ratios.is_empty() {
        let median = ratios[ratios.len() / 2];
        let p10 = ratios[((ratios.len() as f64 * 0.1) as usize).min(ratios.len() - 1)];
        let p90 = ratios[((ratios.len() as f64 * 0.9) as usize).min(ratios.len() - 1)];
        println!("ratio p10 {p10:.2}x  median {median:.2}x  p90 {p90:.2}x");
    }
}
