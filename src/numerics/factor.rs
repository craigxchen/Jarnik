use crate::numerics::exact::isqrt_u128;
use serde::Serialize;
use std::collections::HashMap;

#[derive(Clone, Debug, PartialEq, Eq, Serialize)]
pub struct PrimeFactor {
    pub prime: u64,
    pub exponent: u32,
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize)]
pub struct SplitPrime {
    pub prime: u64,
    pub exponent: u32,
    pub a: i128,
    pub b: i128,
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize)]
pub struct FactorizationData {
    pub n: u64,
    pub factors: Vec<PrimeFactor>,
    pub split_primes: Vec<SplitPrime>,
    pub two_exponent: u32,
    pub inert_scale: u128,
    pub is_sum_of_two_squares: bool,
}

pub fn primes_up_to(limit: u64) -> Vec<u64> {
    if limit < 2 {
        return Vec::new();
    }
    let mut is_composite = vec![false; (limit + 1) as usize];
    let mut primes = Vec::new();
    for candidate in 2..=limit {
        if is_composite[candidate as usize] {
            continue;
        }
        primes.push(candidate);
        if candidate <= limit / candidate {
            let mut multiple = candidate * candidate;
            while multiple <= limit {
                is_composite[multiple as usize] = true;
                multiple += candidate;
            }
        }
    }
    primes
}

pub fn factor_u64(n: u64) -> Vec<PrimeFactor> {
    let limit = isqrt_u128(n as u128) as u64 + 1;
    let primes = primes_up_to(limit);
    factor_u64_with_primes(n, &primes)
}

pub fn factor_u64_with_primes(mut n: u64, primes: &[u64]) -> Vec<PrimeFactor> {
    let mut factors = Vec::new();
    if n < 2 {
        return factors;
    }
    for &prime in primes {
        if prime > n / prime {
            break;
        }
        if n % prime != 0 {
            continue;
        }
        let mut exponent = 0_u32;
        while n % prime == 0 {
            n /= prime;
            exponent += 1;
        }
        factors.push(PrimeFactor { prime, exponent });
    }
    if n > 1 {
        factors.push(PrimeFactor {
            prime: n,
            exponent: 1,
        });
    }
    factors
}

pub fn analyze_factorization(n: u64, primes: &[u64]) -> FactorizationData {
    let mut representation_cache = HashMap::new();
    analyze_factorization_cached(n, primes, &mut representation_cache)
}

pub fn analyze_factorization_cached(
    n: u64,
    primes: &[u64],
    representation_cache: &mut HashMap<u64, Option<(u64, u64)>>,
) -> FactorizationData {
    let factors = factor_u64_with_primes(n, primes);
    let mut split_primes = Vec::new();
    let mut two_exponent = 0_u32;
    let mut inert_scale = 1_u128;
    let mut is_sum_of_two_squares = true;

    for factor in &factors {
        match factor.prime % 4 {
            0 => unreachable!("only 2 is even and 2 % 4 != 0"),
            1 => {
                let representation = *representation_cache
                    .entry(factor.prime)
                    .or_insert_with(|| represent_prime_as_sum_of_squares(factor.prime));
                if let Some((a, b)) = representation {
                    split_primes.push(SplitPrime {
                        prime: factor.prime,
                        exponent: factor.exponent,
                        a: a as i128,
                        b: b as i128,
                    });
                } else {
                    is_sum_of_two_squares = false;
                }
            }
            2 => {
                two_exponent = factor.exponent;
            }
            3 => {
                if factor.exponent % 2 == 1 {
                    is_sum_of_two_squares = false;
                } else {
                    inert_scale *= pow_u128(factor.prime as u128, factor.exponent / 2);
                }
            }
            _ => unreachable!(),
        }
    }

    FactorizationData {
        n,
        factors,
        split_primes,
        two_exponent,
        inert_scale,
        is_sum_of_two_squares,
    }
}

pub fn represent_prime_as_sum_of_squares(p: u64) -> Option<(u64, u64)> {
    if p == 2 {
        return Some((1, 1));
    }
    if p % 4 != 1 {
        return None;
    }
    if let Some(representation) = cornacchia_prime_sum_of_squares(p) {
        return Some(representation);
    }

    // Fallback keeps the routine simple to audit if Cornacchia ever sees an
    // unexpected input. In normal factorization calls, p is prime and this path
    // should not be used.
    let root = isqrt_u128(p as u128) as u64;
    for a in 1..=root {
        let a2 = a * a;
        if a2 > p {
            break;
        }
        let b2 = p - a2;
        let b = isqrt_u128(b2 as u128) as u64;
        if b * b == b2 {
            return if a >= b { Some((a, b)) } else { Some((b, a)) };
        }
    }
    None
}

fn cornacchia_prime_sum_of_squares(p: u64) -> Option<(u64, u64)> {
    let sqrt_minus_one = sqrt_minus_one_mod_prime(p)?;
    let mut r0 = p;
    let mut r1 = sqrt_minus_one.min(p - sqrt_minus_one);
    while (r1 as u128) * (r1 as u128) > p as u128 {
        let r2 = r0 % r1;
        r0 = r1;
        r1 = r2;
    }
    let a = r1;
    let a2 = (a as u128) * (a as u128);
    let b2 = p as u128 - a2;
    let b = isqrt_u128(b2) as u64;
    if (b as u128) * (b as u128) == b2 {
        if a >= b {
            Some((a, b))
        } else {
            Some((b, a))
        }
    } else {
        None
    }
}

fn sqrt_minus_one_mod_prime(p: u64) -> Option<u64> {
    for base in 2..p {
        if mod_pow(base, (p - 1) / 2, p) == p - 1 {
            let root = mod_pow(base, (p - 1) / 4, p);
            if ((root as u128) * (root as u128)) % p as u128 == p as u128 - 1 {
                return Some(root);
            }
        }
    }
    None
}

fn mod_pow(base: u64, mut exponent: u64, modulus: u64) -> u64 {
    let mut acc = 1_u128;
    let modulus_u128 = modulus as u128;
    let mut base_u128 = (base % modulus) as u128;
    while exponent > 0 {
        if exponent & 1 == 1 {
            acc = (acc * base_u128) % modulus_u128;
        }
        exponent >>= 1;
        if exponent > 0 {
            base_u128 = (base_u128 * base_u128) % modulus_u128;
        }
    }
    acc as u64
}

pub fn pow_u128(mut base: u128, mut exponent: u32) -> u128 {
    let mut acc = 1_u128;
    while exponent > 0 {
        if exponent & 1 == 1 {
            acc *= base;
        }
        exponent >>= 1;
        if exponent > 0 {
            base *= base;
        }
    }
    acc
}

#[cfg(test)]
mod tests {
    use super::{
        analyze_factorization, factor_u64, primes_up_to, represent_prime_as_sum_of_squares,
    };

    #[test]
    fn trial_factorization() {
        assert_eq!(
            factor_u64(2 * 2 * 5 * 13 * 13),
            vec![
                super::PrimeFactor {
                    prime: 2,
                    exponent: 2
                },
                super::PrimeFactor {
                    prime: 5,
                    exponent: 1
                },
                super::PrimeFactor {
                    prime: 13,
                    exponent: 2
                }
            ]
        );
    }

    #[test]
    fn split_prime_representatives_are_canonical() {
        assert_eq!(represent_prime_as_sum_of_squares(5), Some((2, 1)));
        assert_eq!(represent_prime_as_sum_of_squares(13), Some((3, 2)));
        assert_eq!(represent_prime_as_sum_of_squares(3), None);
    }

    #[test]
    fn sum_of_two_squares_filter() {
        let primes = primes_up_to(100);
        assert!(analyze_factorization(65, &primes).is_sum_of_two_squares);
        assert!(analyze_factorization(45, &primes).is_sum_of_two_squares);
        assert!(!analyze_factorization(21, &primes).is_sum_of_two_squares);
    }
}
