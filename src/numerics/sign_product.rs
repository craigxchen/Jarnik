use std::collections::{HashMap, HashSet};

#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub struct SignFactor {
    pub r: i128,
    pub t: i128,
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct SignCollision {
    pub coeffs: (i128, i128, i128),
    pub masks: Vec<u64>,
    pub products: Vec<(i128, i128)>,
}

pub fn sign_coefficient_collisions(
    a: i128,
    factors: &[SignFactor],
    min_distinct_vectors: usize,
) -> Vec<SignCollision> {
    if factors.len() > 63 {
        return Vec::new();
    }

    let mut buckets: HashMap<
        (i128, i128, i128),
        (Vec<u64>, Vec<(i128, i128)>, HashSet<(i128, i128)>),
    > = HashMap::new();
    let mask_count = 1_u64 << factors.len();

    for mask in 0..mask_count {
        let coeffs = coefficient_vector(mask, factors);
        let product = gaussian_product(a, mask, factors);
        let entry = buckets
            .entry(coeffs)
            .or_insert_with(|| (Vec::new(), Vec::new(), HashSet::new()));
        if entry.2.insert(product) {
            entry.0.push(mask);
            entry.1.push(product);
        }
    }

    buckets
        .into_iter()
        .filter_map(|(coeffs, (masks, products, _))| {
            if masks.len() >= min_distinct_vectors {
                Some(SignCollision {
                    coeffs,
                    masks,
                    products,
                })
            } else {
                None
            }
        })
        .collect()
}

pub fn search_sign_product_counterexamples(
    a: i128,
    n: usize,
    max_r: i128,
    max_t: i128,
    min_distinct_vectors: usize,
    max_results: usize,
) -> Vec<(Vec<SignFactor>, SignCollision)> {
    let mut candidates = Vec::new();
    for r in 0..=max_r {
        for t in 1..=max_t {
            candidates.push(SignFactor { r, t });
        }
    }

    let mut results = Vec::new();
    let mut chosen = Vec::new();
    search_factor_sets(
        a,
        n,
        min_distinct_vectors,
        max_results,
        &candidates,
        0,
        &mut chosen,
        &mut results,
    );
    results
}

fn search_factor_sets(
    a: i128,
    n: usize,
    min_distinct_vectors: usize,
    max_results: usize,
    candidates: &[SignFactor],
    start: usize,
    chosen: &mut Vec<SignFactor>,
    results: &mut Vec<(Vec<SignFactor>, SignCollision)>,
) {
    if results.len() >= max_results {
        return;
    }
    if chosen.len() == n {
        for collision in sign_coefficient_collisions(a, chosen, min_distinct_vectors) {
            results.push((chosen.clone(), collision));
            if results.len() >= max_results {
                return;
            }
        }
        return;
    }

    for index in start..candidates.len() {
        chosen.push(candidates[index]);
        search_factor_sets(
            a,
            n,
            min_distinct_vectors,
            max_results,
            candidates,
            index + 1,
            chosen,
            results,
        );
        chosen.pop();
        if results.len() >= max_results {
            return;
        }
    }
}

fn coefficient_vector(mask: u64, factors: &[SignFactor]) -> (i128, i128, i128) {
    let mut c0 = 0_i128;
    let mut c1 = 0_i128;
    let mut c2 = 0_i128;
    for (index, factor) in factors.iter().enumerate() {
        let sigma = if (mask >> index) & 1 == 1 { 1 } else { -1 };
        c0 += sigma * factor.t;
        c1 += sigma * factor.t * factor.r;
        c2 += sigma * (3 * factor.t * factor.r * factor.r - factor.t * factor.t * factor.t);
    }
    (c0, c1, c2)
}

fn gaussian_product(a: i128, mask: u64, factors: &[SignFactor]) -> (i128, i128) {
    let mut product = (1_i128, 0_i128);
    for (index, factor) in factors.iter().enumerate() {
        let sigma = if (mask >> index) & 1 == 1 { 1 } else { -1 };
        let z = (a + factor.r, sigma * factor.t);
        product = (
            product.0 * z.0 - product.1 * z.1,
            product.0 * z.1 + product.1 * z.0,
        );
    }
    product
}

#[cfg(test)]
mod tests {
    use super::{sign_coefficient_collisions, SignFactor};

    #[test]
    fn collision_search_detects_balanced_duplicate_coefficients() {
        let factors = [
            SignFactor { r: 0, t: 1 },
            SignFactor { r: 1, t: 2 },
            SignFactor { r: 3, t: 2 },
            SignFactor { r: 4, t: 1 },
        ];
        let collisions = sign_coefficient_collisions(10, &factors, 2);
        assert!(!collisions.is_empty());
    }
}
