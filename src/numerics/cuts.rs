use crate::numerics::factor::SplitPrime;
use crate::numerics::points::LatticePointRecord;
use serde::Serialize;

#[derive(Clone, Debug)]
pub struct CutInstance {
    pub prime: u64,
    pub threshold: u32,
    pub weight: f64,
    pub pattern: Vec<u8>,
    pub variation: usize,
}

impl CutInstance {
    pub fn crosses_gap(&self, gap: usize) -> bool {
        self.pattern[gap] != self.pattern[gap + 1]
    }

    pub fn transition_sign(&self, gap: usize) -> Option<i8> {
        let from = self.pattern[gap] as i8;
        let to = self.pattern[gap + 1] as i8;
        match to - from {
            1 => Some(1),
            -1 => Some(-1),
            _ => None,
        }
    }
}

#[derive(Clone, Debug, Serialize)]
pub struct TopCutRecord {
    pub prime: u64,
    pub threshold: u32,
    pub weight: f64,
    pub pattern: Vec<u8>,
    pub variation: usize,
}

#[derive(Clone, Debug, Serialize)]
pub struct CutDiagnostics {
    pub total_cut_weight: f64,
    pub crossing_cut_weight: f64,
    pub laminar_cut_weight: f64,
    pub crossing_fraction: f64,
    pub variation_mass: f64,
    pub adjacent_cut_distances: Vec<f64>,
    pub top_crossing_cuts: Vec<TopCutRecord>,
}

pub fn enumerate_cuts(
    split_primes: &[SplitPrime],
    cluster_points: &[&LatticePointRecord],
) -> Vec<CutInstance> {
    let mut cuts = Vec::new();
    for (prime_index, split) in split_primes.iter().enumerate() {
        let weight = (split.prime as f64).ln();
        for threshold in 0..split.exponent {
            let pattern = cluster_points
                .iter()
                .map(|point| {
                    if point.exponents[prime_index] > threshold {
                        1
                    } else {
                        0
                    }
                })
                .collect::<Vec<_>>();
            let variation = pattern.windows(2).filter(|pair| pair[0] != pair[1]).count();
            cuts.push(CutInstance {
                prime: split.prime,
                threshold,
                weight,
                pattern,
                variation,
            });
        }
    }
    cuts
}

pub fn analyze_cuts(cuts: &[CutInstance], cluster_len: usize) -> CutDiagnostics {
    let mut total_cut_weight = 0.0;
    let mut crossing_cut_weight = 0.0;
    let mut laminar_cut_weight = 0.0;
    let mut variation_mass = 0.0;
    let mut adjacent_cut_distances = vec![0.0; cluster_len.saturating_sub(1)];

    for cut in cuts {
        total_cut_weight += cut.weight;
        variation_mass += cut.variation as f64 * cut.weight;
        if cut.variation >= 2 {
            crossing_cut_weight += cut.weight;
        } else {
            laminar_cut_weight += cut.weight;
        }
        for (gap, distance) in adjacent_cut_distances.iter_mut().enumerate() {
            if cut.crosses_gap(gap) {
                *distance += cut.weight;
            }
        }
    }

    let mut top_crossing_cuts = cuts
        .iter()
        .filter(|cut| cut.variation >= 2)
        .map(|cut| TopCutRecord {
            prime: cut.prime,
            threshold: cut.threshold,
            weight: cut.weight,
            pattern: cut.pattern.clone(),
            variation: cut.variation,
        })
        .collect::<Vec<_>>();
    top_crossing_cuts.sort_by(|a, b| {
        b.weight
            .partial_cmp(&a.weight)
            .unwrap()
            .then_with(|| b.variation.cmp(&a.variation))
            .then_with(|| a.prime.cmp(&b.prime))
            .then_with(|| a.threshold.cmp(&b.threshold))
    });
    top_crossing_cuts.truncate(16);

    CutDiagnostics {
        total_cut_weight,
        crossing_cut_weight,
        laminar_cut_weight,
        crossing_fraction: if total_cut_weight > 0.0 {
            crossing_cut_weight / total_cut_weight
        } else {
            0.0
        },
        variation_mass,
        adjacent_cut_distances,
        top_crossing_cuts,
    }
}

#[cfg(test)]
mod tests {
    use crate::numerics::factor::{analyze_factorization, primes_up_to};
    use crate::numerics::points::generate_lattice_points_from_factorization;

    use super::{analyze_cuts, enumerate_cuts};

    #[test]
    fn cuts_detect_variation() {
        let primes = primes_up_to(10);
        let factorization = analyze_factorization(25, &primes);
        let points = generate_lattice_points_from_factorization(&factorization);
        let chosen = [(3, 4), (5, 0), (3, -4)]
            .iter()
            .map(|&(x, y)| {
                points
                    .iter()
                    .find(|point| point.x == x && point.y == y)
                    .unwrap()
            })
            .collect::<Vec<_>>();
        let cuts = enumerate_cuts(&factorization.split_primes, &chosen);
        let diagnostics = analyze_cuts(&cuts, chosen.len());
        assert!(diagnostics.total_cut_weight > 0.0);
        assert!(diagnostics.variation_mass > 0.0);
    }
}
