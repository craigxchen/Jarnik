use crate::numerics::exact::isqrt_u128;
use crate::numerics::factor::FactorizationData;
use crate::numerics::gaussian::{units, Gaussian};
use serde::Serialize;
use std::collections::HashMap;

#[derive(Clone, Debug, Serialize)]
pub struct LatticePointRecord {
    pub x: i128,
    pub y: i128,
    pub angle: f64,
    pub exponents: Vec<u32>,
}

impl LatticePointRecord {
    pub fn point_array(&self) -> [i128; 2] {
        [self.x, self.y]
    }
}

pub fn generate_lattice_points_from_factorization(
    factorization: &FactorizationData,
) -> Vec<LatticePointRecord> {
    if factorization.n == 0 || !factorization.is_sum_of_two_squares {
        return Vec::new();
    }

    let mut base = Gaussian::new(factorization.inert_scale as i128, 0);
    if factorization.two_exponent > 0 {
        base = base.mul(Gaussian::new(1, 1).pow(factorization.two_exponent));
    }

    let mut partials = vec![(base, Vec::<u32>::new())];
    for split in &factorization.split_primes {
        let pi = Gaussian::new(split.a, split.b);
        let pi_conj = pi.conjugate();
        let mut next = Vec::new();
        for (partial, exponents) in &partials {
            for exponent in 0..=split.exponent {
                let factor = pi.pow(exponent).mul(pi_conj.pow(split.exponent - exponent));
                let mut next_exponents = exponents.clone();
                next_exponents.push(exponent);
                next.push((partial.mul(factor), next_exponents));
            }
        }
        partials = next;
    }

    let mut by_point: HashMap<(i128, i128), Vec<u32>> = HashMap::new();
    for (partial, exponents) in partials {
        for unit in units() {
            let z = unit.mul(partial);
            by_point
                .entry((z.re, z.im))
                .or_insert_with(|| exponents.clone());
        }
    }

    let mut points: Vec<LatticePointRecord> = by_point
        .into_iter()
        .map(|((x, y), exponents)| LatticePointRecord {
            x,
            y,
            angle: normalized_angle(x, y),
            exponents,
        })
        .collect();
    sort_points_by_angle(&mut points);
    points
}

pub fn direct_lattice_points(n: u64) -> Vec<LatticePointRecord> {
    if n == 0 {
        return Vec::new();
    }
    let radius_floor = isqrt_u128(n as u128) as i128;
    let mut points = Vec::new();
    for x in -radius_floor..=radius_floor {
        let x2 = (x * x) as u128;
        if x2 > n as u128 {
            continue;
        }
        let y2 = n as u128 - x2;
        let y = isqrt_u128(y2);
        if y * y != y2 {
            continue;
        }
        points.push(LatticePointRecord {
            x,
            y: y as i128,
            angle: normalized_angle(x, y as i128),
            exponents: Vec::new(),
        });
        if y != 0 {
            points.push(LatticePointRecord {
                x,
                y: -(y as i128),
                angle: normalized_angle(x, -(y as i128)),
                exponents: Vec::new(),
            });
        }
    }
    sort_points_by_angle(&mut points);
    points
}

pub fn sort_points_by_angle(points: &mut [LatticePointRecord]) {
    points.sort_by(|a, b| {
        a.angle
            .partial_cmp(&b.angle)
            .unwrap()
            .then_with(|| a.x.cmp(&b.x))
            .then_with(|| a.y.cmp(&b.y))
    });
}

fn normalized_angle(x: i128, y: i128) -> f64 {
    let mut angle = (y as f64).atan2(x as f64);
    if angle < 0.0 {
        angle += std::f64::consts::TAU;
    }
    angle
}

#[cfg(test)]
mod tests {
    use crate::numerics::factor::{analyze_factorization, primes_up_to};

    use super::generate_lattice_points_from_factorization;

    #[test]
    fn generates_all_points_on_twenty_five() {
        let primes = primes_up_to(10);
        let factorization = analyze_factorization(25, &primes);
        let points = generate_lattice_points_from_factorization(&factorization);
        let mut coords = points
            .iter()
            .map(|point| (point.x, point.y))
            .collect::<Vec<_>>();
        coords.sort();
        assert_eq!(coords.len(), 12);
        for point in [(3, 4), (4, 3), (5, 0), (0, 5), (-3, -4), (-4, -3)] {
            assert!(coords.contains(&point));
        }
    }

    #[test]
    fn records_split_prime_exponents_for_five() {
        let primes = primes_up_to(10);
        let factorization = analyze_factorization(5, &primes);
        let points = generate_lattice_points_from_factorization(&factorization);
        let point = points
            .iter()
            .find(|point| point.x == 2 && point.y == 1)
            .expect("2+i should be generated from the canonical prime");
        assert_eq!(point.exponents, vec![1]);
    }
}
