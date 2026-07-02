use crate::numerics::certificate::Certificate;
use crate::numerics::exact::{
    div_exact, is_square_i128, mod_i128, mod_inverse_i128, square_root_i128,
};

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct ExtensionCoefficients {
    pub quad_a: i128,
    pub b1: i128,
    pub b0: i128,
    pub a2: i128,
    pub a1: i128,
    pub a0: i128,
}

impl ExtensionCoefficients {
    pub fn from_certificate(cert: &Certificate) -> Self {
        let quad_a = 2 * cert.h * cert.y * cert.z * cert.v;
        let b1 = cert.big_n * cert.v * cert.v - 2 * cert.h * cert.y * cert.z * cert.q;
        let b0 = -2 * cert.h * cert.y * cert.z * cert.z * cert.v;
        let a2 = b1 * b1 - 8 * cert.h * cert.big_n * cert.y * cert.y * cert.z * cert.v * cert.v;
        let a1 = 2 * b1 * b0;
        let a0 = b0 * b0;
        ExtensionCoefficients {
            quad_a,
            b1,
            b0,
            a2,
            a1,
            a0,
        }
    }

    pub fn with_geom_leading(mut self, cert: &Certificate) -> Self {
        let simplified = simplified_leading_coefficient(cert);
        debug_assert_eq!(
            self.a2, simplified,
            "geometric certificate leading discriminant coefficient mismatch"
        );
        self.a2 = simplified;
        self
    }

    pub fn delta(&self, w_det: i128) -> i128 {
        self.a2 * w_det * w_det + self.a1 * w_det + self.a0
    }

    pub fn b_value(&self, w_det: i128) -> i128 {
        self.b1 * w_det + self.b0
    }
}

pub fn simplified_leading_coefficient(cert: &Certificate) -> i128 {
    -4 * cert.y * cert.y * cert.z * cert.z * cert.v * cert.v
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct ExtensionWitness {
    pub w_det: i128,
    pub t: i128,
    pub wp: i128,
    pub cp: i128,
    pub jp: i128,
    pub sqrt_delta: i128,
}

impl ExtensionWitness {
    pub fn to_json_line(&self) -> String {
        format!(
            "{{\"W\":{},\"t\":{},\"wp\":{},\"cp\":{},\"Jp\":{},\"sqrt_delta\":{}}}",
            self.w_det, self.t, self.wp, self.cp, self.jp, self.sqrt_delta
        )
    }
}

#[derive(Clone, Debug)]
pub struct ExtensionSearchConfig {
    pub max_w: Option<i128>,
    pub primes: Vec<u64>,
    pub max_residues: usize,
    pub use_geom_leading: bool,
}

impl Default for ExtensionSearchConfig {
    fn default() -> Self {
        ExtensionSearchConfig {
            max_w: None,
            primes: vec![3, 5, 7, 11, 13, 17],
            max_residues: 1_000_000,
            use_geom_leading: true,
        }
    }
}

#[derive(Clone, Debug)]
pub struct ExtensionSearchReport {
    pub upper_w: i128,
    pub checked_w: u128,
    pub sieve: Option<QuadraticSieve>,
    pub witnesses: Vec<ExtensionWitness>,
}

pub fn find_extensions(
    cert: &Certificate,
    config: &ExtensionSearchConfig,
) -> Result<ExtensionSearchReport, String> {
    let mut coeffs = ExtensionCoefficients::from_certificate(cert);
    if config.use_geom_leading {
        coeffs = coeffs.with_geom_leading(cert);
    }

    let upper_w = match config.max_w {
        Some(max_w) => max_w,
        None => positive_discriminant_upper(&coeffs).ok_or_else(|| {
            "discriminant window is not finite; pass an explicit max_w".to_string()
        })?,
    };

    let sieve = if config.primes.is_empty() {
        None
    } else {
        Some(QuadraticSieve::build(
            coeffs.a2,
            coeffs.a1,
            coeffs.a0,
            &config.primes,
            config.max_residues,
        )?)
    };

    let mut checked_w = 0_u128;
    let mut witnesses = Vec::new();

    if let Some(sieve) = &sieve {
        for residue in &sieve.residues {
            let mut w_det = *residue as i128;
            if w_det == 0 {
                w_det = sieve.modulus as i128;
            }
            while w_det <= upper_w {
                checked_w += 1;
                push_witnesses_for_w(cert, &coeffs, w_det, &mut witnesses);
                w_det += sieve.modulus as i128;
            }
        }
    } else {
        for w_det in 1..=upper_w {
            checked_w += 1;
            push_witnesses_for_w(cert, &coeffs, w_det, &mut witnesses);
        }
    }

    Ok(ExtensionSearchReport {
        upper_w,
        checked_w,
        sieve,
        witnesses,
    })
}

fn push_witnesses_for_w(
    cert: &Certificate,
    coeffs: &ExtensionCoefficients,
    w_det: i128,
    witnesses: &mut Vec<ExtensionWitness>,
) {
    let delta = coeffs.delta(w_det);
    if delta < 0 || !is_square_i128(delta) {
        return;
    }

    let sqrt_delta = square_root_i128(delta).expect("square checked above");
    let b_value = coeffs.b_value(w_det);
    let denom = 2 * coeffs.quad_a;

    for numerator in [-b_value + sqrt_delta, -b_value - sqrt_delta] {
        let t = match div_exact(numerator, denom) {
            Some(t) => t,
            None => continue,
        };
        if t <= cert.z || 4 * t * t > cert.big_n {
            continue;
        }

        let wp = match div_exact(t * cert.v + cert.y * w_det, cert.z) {
            Some(wp) if wp > 0 => wp,
            _ => continue,
        };
        let cp = match div_exact(cert.q * w_det - (t - cert.z) * cert.v, cert.z) {
            Some(cp) if cp > 0 => cp,
            _ => continue,
        };
        let jp = match div_exact(cert.v * w_det * wp, cert.h) {
            Some(jp) if jp > 0 => jp,
            _ => continue,
        };
        if cert.big_n * jp != 2 * cert.y * cert.z * t * cp {
            continue;
        }
        if (t - cert.z) * cert.v * wp <= cert.y * cp * w_det {
            continue;
        }

        witnesses.push(ExtensionWitness {
            w_det,
            t,
            wp,
            cp,
            jp,
            sqrt_delta,
        });
    }
}

pub fn positive_discriminant_upper(coeffs: &ExtensionCoefficients) -> Option<i128> {
    if coeffs.a2 >= 0 {
        return None;
    }

    let mut hi = 1_i128;
    while coeffs.delta(hi) >= 0 {
        if hi > i128::MAX / 4 {
            return None;
        }
        hi *= 2;
    }

    let mut lo = 0_i128;
    while lo + 1 < hi {
        let mid = lo + (hi - lo) / 2;
        if coeffs.delta(mid) >= 0 {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    Some(lo)
}

#[derive(Clone, Debug)]
pub struct QuadraticSieve {
    pub modulus: u128,
    pub residues: Vec<u128>,
    pub prime_count: usize,
}

impl QuadraticSieve {
    pub fn selectivity(&self) -> f64 {
        if self.modulus == 0 {
            0.0
        } else {
            self.residues.len() as f64 / self.modulus as f64
        }
    }

    pub fn build(
        a2: i128,
        a1: i128,
        a0: i128,
        primes: &[u64],
        max_residues: usize,
    ) -> Result<Self, String> {
        let mut modulus = 1_u128;
        let mut residues = vec![0_u128];

        for &prime in primes {
            if prime < 2 {
                return Err(format!("{} is not a prime modulus", prime));
            }
            let prime_u = prime as u128;
            let allowed = allowed_residues_for_prime(a2, a1, a0, prime_u);
            if allowed.is_empty() {
                return Ok(QuadraticSieve {
                    modulus: modulus * prime_u,
                    residues: Vec::new(),
                    prime_count: primes.len(),
                });
            }

            let inv_modulus = mod_inverse_i128((modulus % prime_u) as i128, prime as i128)
                .ok_or_else(|| format!("moduli are not coprime at prime {}", prime))?;
            let mut next = Vec::with_capacity(residues.len() * allowed.len());
            for &residue in &residues {
                for &allowed_residue in &allowed {
                    let difference = mod_i128(
                        allowed_residue as i128 - (residue % prime_u) as i128,
                        prime_u,
                    );
                    let k = (difference * inv_modulus as u128) % prime_u;
                    next.push(residue + modulus * k);
                }
            }

            modulus *= prime_u;
            residues = next;
            residues.sort_unstable();
            residues.dedup();

            if residues.len() > max_residues {
                return Err(format!(
                    "CRT sieve has {} residues modulo {}; raise max_residues or use fewer primes",
                    residues.len(),
                    modulus
                ));
            }
        }

        Ok(QuadraticSieve {
            modulus,
            residues,
            prime_count: primes.len(),
        })
    }
}

pub fn allowed_residues_for_prime(a2: i128, a1: i128, a0: i128, prime: u128) -> Vec<u128> {
    let mut quadratic_residue = vec![false; prime as usize];
    for x in 0..prime {
        quadratic_residue[((x * x) % prime) as usize] = true;
    }

    let mut allowed = Vec::new();
    for w in 0..prime {
        let delta =
            (mod_i128(a2, prime) * w % prime * w + mod_i128(a1, prime) * w + mod_i128(a0, prime))
                % prime;
        if quadratic_residue[delta as usize] {
            allowed.push(w);
        }
    }
    allowed
}

#[cfg(test)]
mod tests {
    use super::{
        allowed_residues_for_prime, find_extensions, simplified_leading_coefficient,
        ExtensionCoefficients, ExtensionSearchConfig, QuadraticSieve,
    };
    use crate::numerics::certificate::Certificate;

    #[test]
    fn known_certificate_has_simplified_leading_coefficient() {
        let cert = Certificate::known_first_endpoint_certificate();
        let coeffs = ExtensionCoefficients::from_certificate(&cert);
        assert_eq!(coeffs.a2, simplified_leading_coefficient(&cert));
    }

    #[test]
    fn known_certificate_has_no_extension_in_discriminant_window() {
        let cert = Certificate::known_first_endpoint_certificate();
        let report = find_extensions(&cert, &ExtensionSearchConfig::default()).unwrap();
        assert!(report.upper_w > 0);
        assert!(report.witnesses.is_empty());
        assert!(report.checked_w < report.upper_w as u128);
    }

    #[test]
    fn residue_sieve_is_consistent_mod_small_primes() {
        let cert = Certificate::known_first_endpoint_certificate();
        let coeffs = ExtensionCoefficients::from_certificate(&cert);
        let allowed = allowed_residues_for_prime(coeffs.a2, coeffs.a1, coeffs.a0, 7);
        assert!(!allowed.is_empty());

        let sieve = QuadraticSieve::build(coeffs.a2, coeffs.a1, coeffs.a0, &[3, 5, 7], 1000)
            .expect("small sieve should build");
        for residue in sieve.residues {
            assert!(
                allowed_residues_for_prime(coeffs.a2, coeffs.a1, coeffs.a0, 7)
                    .contains(&(residue % 7))
            );
        }
    }
}
