use crate::numerics::cuts::CutInstance;
use serde::Serialize;
use std::collections::{BTreeMap, BTreeSet};

const DEFAULT_K_VALUES: [usize; 5] = [2, 3, 4, 5, 6];
const MAX_EXACT_SUBSETS: u128 = 50_000;
const MOSTLY_ANTIPODAL_THRESHOLD: f64 = 2.0 / 3.0;

#[derive(Clone, Debug, Serialize)]
pub struct TransitionPatternDiagnostics {
    pub reports: Vec<TransitionKReport>,
    pub skipped: Vec<SkippedTransitionSearch>,
}

#[derive(Clone, Debug, Serialize)]
pub struct SkippedTransitionSearch {
    pub k: usize,
    pub gap_count: usize,
    pub subset_count: u128,
    pub reason: String,
}

#[derive(Clone, Debug, Serialize)]
pub struct TransitionKReport {
    pub k: usize,
    pub best_gap_subset: Vec<usize>,
    pub common_crossing_weight: f64,
    pub common_crossing_fraction_of_total: f64,
    pub transition_pattern_distribution: BTreeMap<String, f64>,
    pub hyperplane: HyperplaneReport,
    pub antipodal_pair: AntipodalPairReport,
    pub support: PatternSupportReport,
}

#[derive(Clone, Debug, Serialize)]
pub struct HyperplaneReport {
    pub max_zero_hyperplane_fraction: f64,
    pub best_c: Vec<i8>,
    pub zero_patterns: Vec<String>,
    pub exceeds_1_2: bool,
    pub exceeds_0_55: bool,
    pub exceeds_0_6: bool,
    pub exceeds_0_67: bool,
    pub exceeds_0_75: bool,
}

#[derive(Clone, Debug, Serialize)]
pub struct AntipodalPairReport {
    pub max_antipodal_pair_fraction: f64,
    pub mostly_antipodal_pair: bool,
    pub mostly_threshold: f64,
    pub best_pair: Vec<String>,
    pub best_pair_weight: f64,
}

#[derive(Clone, Debug, Serialize)]
pub struct PatternSupportReport {
    pub support_size: usize,
    pub cube_size: usize,
    pub support_fraction: f64,
    pub affine_subspace_or_coset: bool,
    pub min_pairwise_hamming: usize,
    pub max_pairwise_hamming: usize,
    pub near_uniform_l1: f64,
}

pub fn analyze_transition_patterns(
    cuts: &[CutInstance],
    cluster_len: usize,
    total_cut_weight: f64,
) -> TransitionPatternDiagnostics {
    let gap_count = cluster_len.saturating_sub(1);
    let mut reports = Vec::new();
    let mut skipped = Vec::new();

    for k in DEFAULT_K_VALUES {
        if k > gap_count {
            continue;
        }
        let subset_count = binomial(gap_count as u128, k as u128);
        if subset_count > MAX_EXACT_SUBSETS {
            skipped.push(SkippedTransitionSearch {
                k,
                gap_count,
                subset_count,
                reason: format!(
                    "exact scan cap is {} gap subsets; increase the cap in code if needed",
                    MAX_EXACT_SUBSETS
                ),
            });
            continue;
        }

        let mut best_subset = Vec::new();
        let mut best_weight = -1.0_f64;
        let mut best_distribution = BTreeMap::new();
        for subset in combinations(gap_count, k) {
            let distribution = distribution_for_subset(cuts, &subset);
            let weight = distribution.values().sum::<f64>();
            if weight > best_weight {
                best_weight = weight;
                best_subset = subset;
                best_distribution = distribution;
            }
        }

        let common_crossing_weight = best_weight.max(0.0);
        reports.push(TransitionKReport {
            k,
            best_gap_subset: best_subset,
            common_crossing_weight,
            common_crossing_fraction_of_total: if total_cut_weight > 0.0 {
                common_crossing_weight / total_cut_weight
            } else {
                0.0
            },
            hyperplane: analyze_hyperplanes(k, &best_distribution),
            antipodal_pair: analyze_antipodal_pair(&best_distribution),
            support: analyze_support(k, &best_distribution),
            transition_pattern_distribution: best_distribution,
        });
    }

    TransitionPatternDiagnostics { reports, skipped }
}

fn distribution_for_subset(cuts: &[CutInstance], subset: &[usize]) -> BTreeMap<String, f64> {
    let mut distribution = BTreeMap::new();
    'cuts: for cut in cuts {
        let mut pattern = String::with_capacity(subset.len());
        for &gap in subset {
            match cut.transition_sign(gap) {
                Some(1) => pattern.push('+'),
                Some(-1) => pattern.push('-'),
                _ => continue 'cuts,
            }
        }
        *distribution.entry(pattern).or_insert(0.0) += cut.weight;
    }
    distribution
}

fn analyze_hyperplanes(k: usize, distribution: &BTreeMap<String, f64>) -> HyperplaneReport {
    if distribution.is_empty() {
        return HyperplaneReport {
            max_zero_hyperplane_fraction: 0.0,
            best_c: Vec::new(),
            zero_patterns: Vec::new(),
            exceeds_1_2: false,
            exceeds_0_55: false,
            exceeds_0_6: false,
            exceeds_0_67: false,
            exceeds_0_75: false,
        };
    }

    let total = distribution.values().sum::<f64>();
    let mut best_fraction = -1.0_f64;
    let mut best_c = Vec::new();
    let mut best_zero_patterns = Vec::new();
    for coeffs in coefficient_vectors(k) {
        let normalized_coeffs = normalize_coeffs(&coeffs);
        let mut zero_weight = 0.0;
        let mut zero_patterns = Vec::new();
        for (pattern, weight) in distribution {
            if dot_pattern(&normalized_coeffs, pattern) == 0 {
                zero_weight += *weight;
                zero_patterns.push(pattern.clone());
            }
        }
        let fraction = if total > 0.0 {
            zero_weight / total
        } else {
            0.0
        };
        if fraction > best_fraction + 1e-14
            || ((fraction - best_fraction).abs() <= 1e-14
                && prefer_coefficients(&normalized_coeffs, &best_c))
        {
            best_fraction = fraction;
            best_c = normalized_coeffs;
            best_zero_patterns = zero_patterns;
        }
    }

    HyperplaneReport {
        max_zero_hyperplane_fraction: best_fraction.max(0.0),
        best_c,
        zero_patterns: best_zero_patterns,
        exceeds_1_2: best_fraction > 0.5,
        exceeds_0_55: best_fraction > 0.55,
        exceeds_0_6: best_fraction > 0.6,
        exceeds_0_67: best_fraction > 0.67,
        exceeds_0_75: best_fraction > 0.75,
    }
}

fn analyze_antipodal_pair(distribution: &BTreeMap<String, f64>) -> AntipodalPairReport {
    if distribution.is_empty() {
        return AntipodalPairReport {
            max_antipodal_pair_fraction: 0.0,
            mostly_antipodal_pair: false,
            mostly_threshold: MOSTLY_ANTIPODAL_THRESHOLD,
            best_pair: Vec::new(),
            best_pair_weight: 0.0,
        };
    }

    let total = distribution.values().sum::<f64>();
    let mut seen = BTreeSet::new();
    let mut best_pair = Vec::new();
    let mut best_pair_weight = -1.0_f64;
    for (pattern, weight) in distribution {
        let antipode = antipodal_pattern(pattern);
        if seen.contains(pattern) {
            continue;
        }
        seen.insert(pattern.clone());
        seen.insert(antipode.clone());
        let pair_weight = weight + distribution.get(&antipode).copied().unwrap_or(0.0);
        if pair_weight > best_pair_weight {
            best_pair_weight = pair_weight;
            best_pair = if pattern == &antipode {
                vec![pattern.clone()]
            } else {
                vec![pattern.clone(), antipode]
            };
        }
    }

    let fraction = if total > 0.0 {
        best_pair_weight / total
    } else {
        0.0
    };
    AntipodalPairReport {
        max_antipodal_pair_fraction: fraction,
        mostly_antipodal_pair: fraction >= MOSTLY_ANTIPODAL_THRESHOLD,
        mostly_threshold: MOSTLY_ANTIPODAL_THRESHOLD,
        best_pair,
        best_pair_weight: best_pair_weight.max(0.0),
    }
}

fn analyze_support(k: usize, distribution: &BTreeMap<String, f64>) -> PatternSupportReport {
    let support = distribution.keys().cloned().collect::<Vec<_>>();
    let support_size = support.len();
    let cube_size = 1_usize << k;
    let total = distribution.values().sum::<f64>();
    let near_uniform_l1 = if support_size > 0 && total > 0.0 {
        let expected = 1.0 / support_size as f64;
        distribution
            .values()
            .map(|weight| (weight / total - expected).abs())
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

    PatternSupportReport {
        support_size,
        cube_size,
        support_fraction: if cube_size > 0 {
            support_size as f64 / cube_size as f64
        } else {
            0.0
        },
        affine_subspace_or_coset: is_affine_subspace_or_coset(k, &support),
        min_pairwise_hamming: min_hamming,
        max_pairwise_hamming: max_hamming,
        near_uniform_l1,
    }
}

fn combinations(n: usize, k: usize) -> Vec<Vec<usize>> {
    fn rec(
        n: usize,
        k: usize,
        start: usize,
        current: &mut Vec<usize>,
        output: &mut Vec<Vec<usize>>,
    ) {
        if current.len() == k {
            output.push(current.clone());
            return;
        }
        let remaining = k - current.len();
        for value in start..=n - remaining {
            current.push(value);
            rec(n, k, value + 1, current, output);
            current.pop();
        }
    }

    let mut output = Vec::new();
    if k == 0 || k > n {
        return output;
    }
    rec(n, k, 0, &mut Vec::new(), &mut output);
    output
}

fn coefficient_vectors(k: usize) -> Vec<Vec<i8>> {
    fn rec(k: usize, current: &mut Vec<i8>, output: &mut Vec<Vec<i8>>) {
        if current.len() == k {
            if current.iter().any(|&value| value != 0) {
                output.push(current.clone());
            }
            return;
        }
        for value in -2..=2 {
            current.push(value);
            rec(k, current, output);
            current.pop();
        }
    }

    let mut output = Vec::new();
    rec(k, &mut Vec::new(), &mut output);
    output
}

fn normalize_coeffs(coeffs: &[i8]) -> Vec<i8> {
    let gcd = coeffs
        .iter()
        .filter(|&&coeff| coeff != 0)
        .fold(0_i8, |acc, &coeff| gcd_i8(acc, coeff.abs()))
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

fn prefer_coefficients(candidate: &[i8], incumbent: &[i8]) -> bool {
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
    coeffs: &[i8],
) -> (usize, std::cmp::Reverse<usize>, std::cmp::Reverse<i8>) {
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

fn gcd_i8(mut a: i8, mut b: i8) -> i8 {
    a = a.abs();
    b = b.abs();
    while b != 0 {
        let r = a % b;
        a = b;
        b = r;
    }
    a
}

fn dot_pattern(coeffs: &[i8], pattern: &str) -> i32 {
    coeffs
        .iter()
        .zip(pattern.chars())
        .map(|(&coeff, sign)| {
            let sign_value = if sign == '+' { 1_i32 } else { -1_i32 };
            coeff as i32 * sign_value
        })
        .sum()
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

fn binomial(n: u128, k: u128) -> u128 {
    if k > n {
        return 0;
    }
    let k = k.min(n - k);
    let mut acc = 1_u128;
    for i in 0..k {
        acc = acc * (n - i) / (i + 1);
    }
    acc
}

#[cfg(test)]
mod tests {
    use super::{analyze_antipodal_pair, analyze_hyperplanes, is_affine_subspace_or_coset};
    use std::collections::BTreeMap;

    #[test]
    fn hyperplane_search_finds_balanced_coordinate() {
        let mut distribution = BTreeMap::new();
        distribution.insert("++".to_string(), 1.0);
        distribution.insert("+-".to_string(), 1.0);
        distribution.insert("-+".to_string(), 1.0);
        distribution.insert("--".to_string(), 1.0);
        let report = analyze_hyperplanes(2, &distribution);
        assert_eq!(report.max_zero_hyperplane_fraction, 0.5);
    }

    #[test]
    fn detects_affine_supports() {
        assert!(is_affine_subspace_or_coset(
            3,
            &["+--".to_string(), "+++".to_string()]
        ));
        assert!(!is_affine_subspace_or_coset(
            3,
            &["+++".to_string(), "+--".to_string(), "-+-".to_string()]
        ));
    }

    #[test]
    fn reports_mostly_antipodal_pair() {
        let mut distribution = BTreeMap::new();
        distribution.insert("++-".to_string(), 3.0);
        distribution.insert("--+".to_string(), 2.0);
        distribution.insert("+-+".to_string(), 1.0);
        let report = analyze_antipodal_pair(&distribution);
        assert_eq!(report.best_pair, vec!["++-".to_string(), "--+".to_string()]);
        assert!((report.max_antipodal_pair_fraction - 5.0 / 6.0).abs() < 1e-12);
        assert!(report.mostly_antipodal_pair);
    }

    #[test]
    fn antipodal_pair_handles_one_sided_support() {
        let mut distribution = BTreeMap::new();
        distribution.insert("-+".to_string(), 1.0);
        let report = analyze_antipodal_pair(&distribution);
        assert_eq!(report.max_antipodal_pair_fraction, 1.0);
        assert_eq!(report.best_pair, vec!["-+".to_string(), "+-".to_string()]);
        assert!(report.mostly_antipodal_pair);
    }

    #[test]
    fn hyperplane_tie_break_prefers_endpoint_relation() {
        let mut distribution = BTreeMap::new();
        distribution.insert("+-+".to_string(), 3.0);
        distribution.insert("-+-".to_string(), 2.0);
        let report = analyze_hyperplanes(3, &distribution);
        assert_eq!(report.max_zero_hyperplane_fraction, 1.0);
        assert_eq!(report.best_c, vec![1, 0, -1]);
    }
}
