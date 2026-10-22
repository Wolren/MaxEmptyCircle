//! Randomised parity sweep against GEOS, behind `bench-geos`.

use max_empty_circle::geo_types::{
    Coord, Geometry, LineString, MultiLineString, MultiPoint, Point, Polygon,
};
use max_empty_circle::mec::input::ObstacleSet;
use max_empty_circle::mec::solver::bench_geos::{geos_tolerance_for, solve_geos};
use max_empty_circle::{largest_empty_circle, MecOptions};

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

fn main() {
    if !cfg!(feature = "bench-geos") {
        println!("run with --features bench-geos");
        return;
    }
    let mut worst_pts = (0.0f64, 0u64);
    let mut worst_segs = (0.0f64, 0u64);
    let mut fails = 0usize;
    let mut total = 0usize;

    // Point sets: many seeds, varying count.
    for seed in 1..=400u64 {
        let mut rng = Rng(seed.wrapping_mul(0x9E37_79B9_7F4A_7C15) | 1);
        let n = 3 + (seed as usize % 40);
        let pts: Vec<Point<f64>> = (0..n)
            .map(|_| Point::new(rng.range(0.0, 100.0), rng.range(0.0, 100.0)))
            .collect();
        let g = Geometry::MultiPoint(MultiPoint::new(pts));
        let Ok(set) = ObstacleSet::from_geometry(&g) else {
            continue;
        };
        let Ok(nat) = largest_empty_circle(&g, None, &MecOptions::default()) else {
            continue;
        };
        let Ok(ge) = solve_geos(&set, None) else {
            continue;
        };
        let (a, b, c, d) = set.bounds().unwrap();
        let diag = (c - a).hypot(d - b).max(1.0);
        let tol = (diag * 1e-6).max(geos_tolerance_for(&set) * 10.0);
        let rel = (ge.radius - nat.radius) / ge.radius.max(1e-300);
        total += 1;
        if nat.radius + tol < ge.radius {
            fails += 1;
            println!(
                "PT FAIL seed={seed} n={n} native={:.9} geos={:.9} gap={:.3e}",
                nat.radius,
                ge.radius,
                ge.radius - nat.radius
            );
        }
        if rel > worst_pts.0 {
            worst_pts = (rel, seed);
        }
    }

    // Segment sets.
    for seed in 1..=200u64 {
        let mut rng = Rng(seed.wrapping_mul(0x2545_F491_4F6C_DD1D) | 1);
        let n = 3 + (seed as usize % 20);
        let segs: Vec<LineString<f64>> = (0..n)
            .map(|_| {
                let x = rng.range(0.0, 100.0);
                let y = rng.range(0.0, 100.0);
                LineString::from(vec![
                    Coord { x, y },
                    Coord {
                        x: x + rng.range(-25.0, 25.0),
                        y: y + rng.range(-25.0, 25.0),
                    },
                ])
            })
            .collect();
        let g = Geometry::MultiLineString(MultiLineString::new(segs));
        let Ok(set) = ObstacleSet::from_geometry(&g) else {
            continue;
        };
        let Ok(nat) = largest_empty_circle(&g, None, &MecOptions::default()) else {
            continue;
        };
        let Ok(ge) = solve_geos(&set, None) else {
            continue;
        };
        let (a, b, c, d) = set.bounds().unwrap();
        let diag = (c - a).hypot(d - b).max(1.0);
        let tol = (diag * 1e-5).max(geos_tolerance_for(&set) * 10.0);
        let rel = (ge.radius - nat.radius) / ge.radius.max(1e-300);
        total += 1;
        if nat.radius + tol < ge.radius {
            fails += 1;
            println!(
                "SEG FAIL seed={seed} n={n} native={:.9} geos={:.9} gap={:.3e}",
                nat.radius,
                ge.radius,
                ge.radius - nat.radius
            );
        }
        if rel > worst_segs.0 {
            worst_segs = (rel, seed);
        }
    }

    // Boundary-confined variants.
    let mut fails_b = 0usize;
    for seed in 1..=100u64 {
        let mut rng = Rng(seed.wrapping_mul(0xDEAD_BEEF_1234_5678) | 1);
        let n = 6 + (seed as usize % 30);
        let pts: Vec<Point<f64>> = (0..n)
            .map(|_| Point::new(rng.range(0.0, 100.0), rng.range(0.0, 100.0)))
            .collect();
        let g = Geometry::MultiPoint(MultiPoint::new(pts));
        let x0 = rng.range(10.0, 40.0);
        let y0 = rng.range(10.0, 40.0);
        let side = rng.range(20.0, 50.0);
        let b = Geometry::Polygon(Polygon::new(
            LineString::from(vec![
                (x0, y0),
                (x0 + side, y0),
                (x0 + side, y0 + side),
                (x0, y0 + side),
                (x0, y0),
            ]),
            vec![],
        ));
        let Ok(set) = ObstacleSet::from_geometry(&g) else {
            continue;
        };
        let Ok(nat) = largest_empty_circle(&g, Some(&b), &MecOptions::default()) else {
            continue;
        };
        let Ok(ge) = solve_geos(&set, Some(&b)) else {
            continue;
        };
        let tol = 0.01;
        total += 1;
        if nat.radius + tol < ge.radius {
            fails_b += 1;
            println!(
                "BND FAIL seed={seed} n={n} native={:.9} geos={:.9}",
                nat.radius, ge.radius
            );
        }
    }

    println!("\n=== SUMMARY");
    println!("cases: {total}");
    println!("fails: {}  (boundary: {fails_b})", fails + fails_b);
    println!(
        "worst point rel gap: {:.3e} (seed {})",
        worst_pts.0, worst_pts.1
    );
    println!(
        "worst seg   rel gap: {:.3e} (seed {})",
        worst_segs.0, worst_segs.1
    );
}
