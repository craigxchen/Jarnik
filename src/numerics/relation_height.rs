use crate::numerics::cuts::CutInstance;
use serde::Serialize;
use std::collections::{BTreeMap, BTreeSet};

const FULL_RELATION_L1_BUDGETS: [usize; 9] = [1, 2, 3, 4, 5, 6, 8, 10, 12];

#[derive(Clone, Debug, Serialize)]
pub struct FullRelationHeightDiagnostics {
    pub height_model: String,
    pub budgets: Vec<FullRelationBudgetReport>,
    pub best_overall: Option<FullRelationBudgetReport>,
    pub balanced_sum_zero_budgets: Vec<FullRelationBudgetReport>,
    pub best_balanced_sum_zero: Option<FullRelationBudgetReport>,
    pub full_transition_shape: FullTransitionShapeReport,
}

#[derive(Clone, Debug, Serialize)]
pub struct FullRelationBudgetReport {
    pub coefficient_l1_budget: usize,
    pub min_h: f64,
    pub min_h_over_log_r: f64,
    pub min_h_over_total_cut_weight: f64,
    pub best_c: Vec<i32>,
    pub best_c_l1: usize,
    pub zero_cuts: Vec<CutTableRecord>,
    pub zero_cut_weight: f64,
    pub zero_cut_fraction: f64,
    pub local_second_difference_like: bool,
    pub local_second_difference_support: Option<[usize; 3]>,
}

#[derive(Clone, Debug, Serialize)]
pub struct CutTableRecord {
    pub prime: u64,
    pub threshold: u32,
    pub weight: f64,
    pub pattern: Vec<u8>,
    pub variation: usize,
}

#[derive(Clone, Debug, Serialize)]
pub struct FullTransitionShapeReport {
    pub gap_count: usize,
    pub total_full_crossing_weight: f64,
    pub full_crossing_fraction_of_total_cut_weight: f64,
    pub pattern_distribution: BTreeMap<String, f64>,
    pub support_size: usize,
    pub cube_size: usize,
    pub support_fraction: f64,
    pub affine_subspace_or_coset: bool,
    pub hadamard_like: bool,
    pub antipodal_pair: bool,
    pub min_pairwise_hamming: usize,
    pub max_pairwise_hamming: usize,
    pub near_uniform_l1: f64,
}

pub fn analyze_full_relation_height(
    cuts: &[CutInstance],
    cluster_len: usize,
    total_cut_weight: f64,
    log_r: f64,
) -> FullRelationHeightDiagnostics {
    let budgets = FULL_RELATION_L1_BUDGETS
        .iter()
        .filter_map(|&budget| {
            analyze_budget(cuts, cluster_len, total_cut_weight, log_r, budget, false)
        })
        .collect::<Vec<_>>();

    let balanced_sum_zero_budgets = FULL_RELATION_L1_BUDGETS
        .iter()
        .filter_map(|&budget| {
            analyze_budget(cuts, cluster_len, total_cut_weight, log_r, budget, true)
        })
        .collect::<Vec<_>>();

    let best_overall = budgets.iter().cloned().min_by(|left, right| {
        left.min_h
            .partial_cmp(&right.min_h)
            .unwrap()
            .then_with(|| left.best_c_l1.cmp(&right.best_c_l1))
            .then_with(|| left.coefficient_l1_budget.cmp(&right.coefficient_l1_budget))
            .then_with(|| left.best_c.cmp(&right.best_c))
    });
    let best_balanced_sum_zero = balanced_sum_zero_budgets
        .iter()
        .cloned()
        .min_by(|left, right| {
            left.min_h
                .partial_cmp(&right.min_h)
                .unwrap()
                .then_with(|| left.best_c_l1.cmp(&right.best_c_l1))
                .then_with(|| left.coefficient_l1_budget.cmp(&right.coefficient_l1_budget))
                .then_with(|| left.best_c.cmp(&right.best_c))
        });

    FullRelationHeightDiagnostics {
        height_model: "canonical_and_conjugate_layer_cuts".to_string(),
        budgets,
        best_overall,
        balanced_sum_zero_budgets,
        best_balanced_sum_zero,
        full_transition_shape: analyze_full_transition_shape(cuts, cluster_len, total_cut_weight),
    }
}

fn analyze_budget(
    cuts: &[CutInstance],
    cluster_len: usize,
    total_cut_weight: f64,
    log_r: f64,
    budget: usize,
    require_zero_sum: bool,
) -> Option<FullRelationBudgetReport> {
    let mut best_h = f64::INFINITY;
    let mut best_c = Vec::new();
    let mut best_zero_cuts = Vec::new();
    let mut best_zero_weight = 0.0;

    for coeffs in coefficient_vectors_l1(cluster_len, budget) {
        if require_zero_sum && coeffs.iter().sum::<i32>() != 0 {
            continue;
        }
        let (h, zero_weight, zero_cuts) =
            relation_height_for_coeffs(&coeffs, cuts, total_cut_weight);
        if h < best_h - 1e-14
            || ((h - best_h).abs() <= 1e-14 && prefer_coefficients(&coeffs, &best_c))
        {
            best_h = h;
            best_c = coeffs;
            best_zero_cuts = zero_cuts;
            best_zero_weight = zero_weight;
        }
    }

    if !best_h.is_finite() {
        return None;
    }

    let min_h = best_h;
    let (local_second_difference_like, local_second_difference_support) =
        second_difference_support(&best_c);
    Some(FullRelationBudgetReport {
        coefficient_l1_budget: budget,
        min_h,
        min_h_over_log_r: if log_r > 0.0 { min_h / log_r } else { 0.0 },
        min_h_over_total_cut_weight: if total_cut_weight > 0.0 {
            min_h / total_cut_weight
        } else {
            0.0
        },
        best_c_l1: coeff_l1(&best_c),
        best_c,
        zero_cuts: best_zero_cuts,
        zero_cut_weight: best_zero_weight,
        zero_cut_fraction: if total_cut_weight > 0.0 {
            best_zero_weight / total_cut_weight
        } else {
            0.0
        },
        local_second_difference_like,
        local_second_difference_support,
    })
}

fn relation_height_for_coeffs(
    coeffs: &[i32],
    cuts: &[CutInstance],
    _total_cut_weight: f64,
) -> (f64, f64, Vec<CutTableRecord>) {
    let mut height = 0.0;
    let mut zero_weight = 0.0;
    let mut zero_cuts = Vec::new();
    let coefficient_sum = coeffs.iter().sum::<i32>();
    for cut in cuts {
        let dot = dot_indicator(coeffs, &cut.pattern);
        let conjugate_dot = coefficient_sum - dot;
        height += (dot.abs() + conjugate_dot.abs()) as f64 * cut.weight;
        if dot == 0 && conjugate_dot == 0 {
            zero_weight += cut.weight;
            zero_cuts.push(CutTableRecord {
                prime: cut.prime,
                threshold: cut.threshold,
                weight: cut.weight,
                pattern: cut.pattern.clone(),
                variation: cut.variation,
            });
        }
    }
    (height, zero_weight, zero_cuts)
}

fn analyze_full_transition_shape(
    cuts: &[CutInstance],
    cluster_len: usize,
    total_cut_weight: f64,
) -> FullTransitionShapeReport {
    let gap_count = cluster_len.saturating_sub(1);
    let mut distribution = BTreeMap::new();
    for cut in cuts {
        if cut.variation != gap_count {
            continue;
        }
        let pattern = (0..gap_count)
            .map(|gap| {
                if cut.pattern[gap + 1] > cut.pattern[gap] {
                    '+'
                } else {
                    '-'
                }
            })
            .collect::<String>();
        *distribution.entry(pattern).or_insert(0.0) += cut.weight;
    }

    let support = distribution.keys().cloned().collect::<Vec<_>>();
    let support_size = support.len();
    let cube_size = if gap_count < usize::BITS as usize {
        1_usize << gap_count
    } else {
        0
    };
    let total_full_crossing_weight = distribution.values().sum::<f64>();
    let near_uniform_l1 = if support_size > 0 && total_full_crossing_weight > 0.0 {
        let expected = 1.0 / support_size as f64;
        distribution
            .values()
            .map(|weight| (weight / total_full_crossing_weight - expected).abs())
            .sum()
    } else {
        0.0
    };

    let mut min_hamming = usize::MAX;
    let mut max_hamming = 0_usize;
    for i in 0..support.len() {
        for j in (i + 1)..support.len() {
            let distance = hamming(&support[i], &support[j]);
            min_hamming = min_hamming.min(distance);
            max_hamming = max_hamming.max(distance);
        }
    }
    if min_hamming == usize::MAX {
        min_hamming = 0;
    }

    let affine_subspace_or_coset = is_affine_subspace_or_coset(gap_count, &support);
    let antipodal_pair = support_size == 2
        && support
            .first()
            .map(|pattern| antipodal_pattern(pattern) == support[1])
            .unwrap_or(false);
    let hadamard_like = affine_subspace_or_coset
        && support_size > 0
        && (antipodal_pair || (min_hamming * 2 >= gap_count && near_uniform_l1 <= 0.25));

    FullTransitionShapeReport {
        gap_count,
        total_full_crossing_weight,
        full_crossing_fraction_of_total_cut_weight: if total_cut_weight > 0.0 {
            total_full_crossing_weight / total_cut_weight
        } else {
            0.0
        },
        pattern_distribution: distribution,
        support_size,
        cube_size,
        support_fraction: if cube_size > 0 {
            support_size as f64 / cube_size as f64
        } else {
            0.0
        },
        affine_subspace_or_coset,
        hadamard_like,
        antipodal_pair,
        min_pairwise_hamming: min_hamming,
        max_pairwise_hamming: max_hamming,
        near_uniform_l1,
    }
}

fn coefficient_vectors_l1(len: usize, budget: usize) -> Vec<Vec<i32>> {
    fn rec(
        len: usize,
        index: usize,
        remaining: usize,
        current: &mut Vec<i32>,
        output: &mut BTreeSet<Vec<i32>>,
    ) {
        if index == len {
            if current.iter().any(|&value| value != 0) {
                output.insert(normalize_coeffs(current));
            }
            return;
        }

        let limit = remaining as i32;
        for value in -limit..=limit {
            let used = value.unsigned_abs() as usize;
            current.push(value);
            rec(len, index + 1, remaining - used, current, output);
            current.pop();
        }
    }

    let mut output = BTreeSet::new();
    rec(len, 0, budget, &mut Vec::new(), &mut output);
    output.into_iter().collect()
}

fn normalize_coeffs(coeffs: &[i32]) -> Vec<i32> {
    let gcd = coeffs
        .iter()
        .filter(|&&coeff| coeff != 0)
        .fold(0_i32, |acc, &coeff| gcd_i32(acc, coeff.abs()))
        .max(1);
    let mut normalized = coeffs.iter().map(|&coeff| coeff / gcd).collect::<Vec<_>>();
    if normalized
        .iter()
        .find(|&&coeff| coeff != 0)
        .map(|&coeff| coeff < 0)
        .unwrap_or(false)
    {
        for coeff in &mut normalized {
            *coeff = -*coeff;
        }
    }
    normalized
}

fn prefer_coefficients(candidate: &[i32], incumbent: &[i32]) -> bool {
    if incumbent.is_empty() {
        return true;
    }
    let candidate_score = coefficient_tie_score(candidate);
    let incumbent_score = coefficient_tie_score(incumbent);
    if candidate_score != incumbent_score {
        return candidate_score > incumbent_score;
    }
    candidate < incumbent
}

fn coefficient_tie_score(
    coeffs: &[i32],
) -> (usize, std::cmp::Reverse<usize>, std::cmp::Reverse<i32>) {
    let nonzero_indices = coeffs
        .iter()
        .enumerate()
        .filter_map(|(idx, &coeff)| if coeff != 0 { Some(idx) } else { None })
        .collect::<Vec<_>>();
    let span = match (nonzero_indices.first(), nonzero_indices.last()) {
        (Some(first), Some(last)) => last - first,
        _ => 0,
    };
    let nonzero_count = nonzero_indices.len();
    let max_abs = coeffs.iter().map(|coeff| coeff.abs()).max().unwrap_or(0);
    (
        span,
        std::cmp::Reverse(nonzero_count),
        std::cmp::Reverse(max_abs),
    )
}

fn second_difference_support(coeffs: &[i32]) -> (bool, Option<[usize; 3]>) {
    for start in 0..coeffs.len().saturating_sub(2) {
        let mut expected = vec![0_i32; coeffs.len()];
        expected[start] = 1;
        expected[start + 1] = -2;
        expected[start + 2] = 1;
        if coeffs == expected {
            return (true, Some([start, start + 1, start + 2]));
        }
    }
    (false, None)
}

fn dot_indicator(coeffs: &[i32], pattern: &[u8]) -> i32 {
    coeffs
        .iter()
        .zip(pattern.iter())
        .map(|(&coeff, &indicator)| coeff * indicator as i32)
        .sum()
}

fn coeff_l1(coeffs: &[i32]) -> usize {
    coeffs
        .iter()
        .map(|coeff| coeff.unsigned_abs() as usize)
        .sum()
}

fn gcd_i32(mut a: i32, mut b: i32) -> i32 {
    a = a.abs();
    b = b.abs();
    while b != 0 {
        let r = a % b;
        a = b;
        b = r;
    }
    a
}

fn hamming(a: &str, b: &str) -> usize {
    a.chars()
        .zip(b.chars())
        .filter(|(left, right)| left != right)
        .count()
}

fn antipodal_pattern(pattern: &str) -> String {
    pattern
        .chars()
        .map(|ch| if ch == '+' { '-' } else { '+' })
        .collect()
}

fn is_affine_subspace_or_coset(k: usize, support: &[String]) -> bool {
    if support.is_empty() || !support.len().is_power_of_two() {
        return false;
    }
    let base = pattern_bits(&support[0]);
    let shifted = support
        .iter()
        .map(|pattern| pattern_bits(pattern) ^ base)
        .collect::<BTreeSet<_>>();
    if shifted.len() != support.len() || !shifted.contains(&0) {
        return false;
    }
    let mask_limit = 1_u128 << k;
    if shifted.iter().any(|&bits| bits >= mask_limit) {
        return false;
    }
    for &a in &shifted {
        for &b in &shifted {
            if !shifted.contains(&(a ^ b)) {
                return false;
            }
        }
    }
    true
}

fn pattern_bits(pattern: &str) -> u128 {
    pattern.chars().enumerate().fold(0_u128, |acc, (idx, ch)| {
        if ch == '-' {
            acc | (1_u128 << idx)
        } else {
            acc
        }
    })
}

#[cfg(test)]
mod tests {
    use super::{analyze_full_relation_height, second_difference_support};
    use crate::numerics::cuts::CutInstance;

    #[test]
    fn detects_local_second_difference() {
        assert_eq!(
            second_difference_support(&[0, 1, -2, 1]),
            (true, Some([1, 2, 3]))
        );
        assert_eq!(second_difference_support(&[1, 0, -1]), (false, None));
    }

    #[test]
    fn full_relation_budget_annihilates_full_cut_patterns() {
        let cuts = vec![
            CutInstance {
                prime: 5,
                threshold: 0,
                weight: 1.0,
                pattern: vec![0, 1, 0],
                variation: 2,
            },
            CutInstance {
                prime: 13,
                threshold: 0,
                weight: 2.0,
                pattern: vec![1, 0, 1],
                variation: 2,
            },
        ];
        let report = analyze_full_relation_height(&cuts, 3, 3.0, 10.0);
        let budget_two = report
            .budgets
            .iter()
            .find(|budget| budget.coefficient_l1_budget == 2)
            .unwrap();
        assert_eq!(budget_two.min_h, 0.0);
        assert_eq!(budget_two.best_c, vec![1, 0, -1]);
    }
}
