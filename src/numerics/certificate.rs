use crate::numerics::exact::{div_exact, gcd_i128, mod_i128, mod_inverse_i128, square_root_i128};

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Certificate {
    pub x: i128,
    pub y: i128,
    pub z: i128,
    pub p: i128,
    pub q: i128,
    pub r: i128,
    pub u: i128,
    pub v: i128,
    pub w: i128,
    pub c: i128,
    pub j: i128,
    pub h: i128,
    pub big_n: i128,
}

impl Certificate {
    pub fn from_root_data(x: i128, y: i128, z: i128, u: i128, v: i128) -> Option<Self> {
        if !(0 < x && x < y && y < z && u > 0 && v > 0) {
            return None;
        }

        let p = y - x;
        let q = z - y;
        let r = z - x;
        let w = div_exact(z * u + x * v, y)?;
        let c = div_exact(p * v - q * u, y)?;
        if w <= 0 || c <= 0 {
            return None;
        }

        let n2_left = p * v * w - z * c * u;
        let n2_right = q * u * w + x * c * v;
        if n2_left != n2_right {
            return None;
        }

        let n3 = q * u * w - x * c * v;
        if n3 <= 0 {
            return None;
        }

        let j_squared = 4 * x * z * c * w * u * v - n2_left * n2_left;
        let j = square_root_i128(j_squared)?;
        if j <= 0 {
            return None;
        }

        let h = div_exact(u * v * w, j)?;
        let big_n = div_exact(2 * x * y * z * c, j)?;
        if h <= 0 || big_n <= 0 {
            return None;
        }

        let cert = Certificate {
            x,
            y,
            z,
            p,
            q,
            r,
            u,
            v,
            w,
            c,
            j,
            h,
            big_n,
        };

        if !cert.is_endpoint_valid() || !cert.reconstructs() || !cert.geom_square_checks() {
            return None;
        }
        Some(cert)
    }

    pub fn known_first_endpoint_certificate() -> Self {
        Certificate::from_root_data(233, 1098, 1597, 6, 6)
            .expect("known certificate should satisfy the certificate equations")
    }

    pub fn n2(&self) -> i128 {
        self.p * self.v * self.w - self.z * self.c * self.u
    }

    pub fn n2_alt(&self) -> i128 {
        self.q * self.u * self.w + self.x * self.c * self.v
    }

    pub fn n3(&self) -> i128 {
        self.q * self.u * self.w - self.x * self.c * self.v
    }

    pub fn is_endpoint_valid(&self) -> bool {
        4 * self.z * self.z <= self.big_n
    }

    pub fn radius_squared(&self) -> i128 {
        self.h * self.h * self.big_n
    }

    pub fn reconstructs(&self) -> bool {
        self.p == self.y - self.x
            && self.q == self.z - self.y
            && self.r == self.z - self.x
            && self.y * self.w == self.z * self.u + self.x * self.v
            && self.y * self.c == self.p * self.v - self.q * self.u
            && self.h * self.j == self.u * self.v * self.w
            && self.big_n * self.j == 2 * self.x * self.y * self.z * self.c
    }

    pub fn geom_square_checks(&self) -> bool {
        let n2 = self.n2();
        let n3 = self.n3();
        self.n2_alt() == n2
            && self.j * self.j + n2 * n2 == 4 * self.x * self.z * self.c * self.w * self.u * self.v
            && self.j * self.j + n3 * n3 == 4 * self.x * self.y * self.c * self.u * self.v * self.w
    }

    pub fn to_json_line(&self) -> String {
        format!(
            "{{\"x\":{},\"y\":{},\"z\":{},\"p\":{},\"q\":{},\"r\":{},\"U\":{},\"V\":{},\"w\":{},\"c\":{},\"J\":{},\"h\":{},\"N\":{}}}",
            self.x,
            self.y,
            self.z,
            self.p,
            self.q,
            self.r,
            self.u,
            self.v,
            self.w,
            self.c,
            self.j,
            self.h,
            self.big_n
        )
    }
}

#[derive(Clone, Debug)]
pub struct CertificateSearch {
    pub min_z: i128,
    pub max_z: i128,
    pub max_det: i128,
    pub limit: Option<usize>,
}

#[derive(Clone, Debug)]
pub struct CertificateSearchReport {
    pub checked_blocks: u128,
    pub congruence_hits: u128,
    pub certificates: Vec<Certificate>,
}

impl CertificateSearch {
    pub fn run(&self) -> CertificateSearchReport {
        let mut certificates = Vec::new();
        let mut checked_blocks = 0_u128;
        let mut congruence_hits = 0_u128;

        'outer: for z in self.min_z.max(3)..=self.max_z {
            for y in 2..z {
                for u in 1..=self.max_det {
                    for v in 1..=self.max_det {
                        checked_blocks += 1;

                        let rhs = -z * u;
                        let gcd = gcd_i128(v, y);
                        if rhs % gcd != 0 {
                            continue;
                        }

                        let reduced_modulus = y / gcd;
                        let residue = if reduced_modulus == 1 {
                            0
                        } else {
                            let reduced_v = v / gcd;
                            let reduced_rhs = rhs / gcd;
                            let inverse = match mod_inverse_i128(
                                mod_i128(reduced_v, reduced_modulus as u128) as i128,
                                reduced_modulus,
                            ) {
                                Some(inverse) => inverse,
                                None => continue,
                            };
                            ((mod_i128(reduced_rhs, reduced_modulus as u128) as i128 * inverse)
                                % reduced_modulus
                                + reduced_modulus)
                                % reduced_modulus
                        };

                        for step in 0..gcd {
                            let x = residue + step * reduced_modulus;
                            if x <= 0 || x >= y {
                                continue;
                            }
                            congruence_hits += 1;
                            if let Some(cert) = Certificate::from_root_data(x, y, z, u, v) {
                                certificates.push(cert);
                                if self
                                    .limit
                                    .map_or(false, |limit| certificates.len() >= limit)
                                {
                                    break 'outer;
                                }
                            }
                        }
                    }
                }
            }
        }

        CertificateSearchReport {
            checked_blocks,
            congruence_hits,
            certificates,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::{Certificate, CertificateSearch};

    #[test]
    fn known_certificate_is_recovered() {
        let cert = Certificate::known_first_endpoint_certificate();
        assert_eq!(cert.x, 233);
        assert_eq!(cert.y, 1098);
        assert_eq!(cert.z, 1597);
        assert_eq!(cert.u, 6);
        assert_eq!(cert.v, 6);
        assert_eq!(cert.w, 10);
        assert_eq!(cert.c, 2);
        assert_eq!(cert.j, 72);
        assert_eq!(cert.h, 5);
        assert_eq!(cert.big_n, 22_698_161);
        assert_eq!(cert.radius_squared(), 567_454_025);
        assert!(cert.reconstructs());
        assert!(cert.geom_square_checks());
        assert!(cert.is_endpoint_valid());
    }

    #[test]
    fn small_search_finds_known_certificate() {
        let report = CertificateSearch {
            min_z: 1597,
            max_z: 1597,
            max_det: 6,
            limit: None,
        }
        .run();
        assert!(report
            .certificates
            .iter()
            .any(|cert| cert == &Certificate::known_first_endpoint_certificate()));
    }
}
