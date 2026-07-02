use serde::Serialize;

#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, Serialize)]
pub struct Gaussian {
    pub re: i128,
    pub im: i128,
}

impl Gaussian {
    pub fn new(re: i128, im: i128) -> Self {
        Self { re, im }
    }

    pub fn zero() -> Self {
        Self { re: 0, im: 0 }
    }

    pub fn one() -> Self {
        Self { re: 1, im: 0 }
    }

    pub fn conjugate(self) -> Self {
        Self {
            re: self.re,
            im: -self.im,
        }
    }

    pub fn norm(self) -> i128 {
        self.re * self.re + self.im * self.im
    }

    pub fn mul(self, rhs: Self) -> Self {
        Self {
            re: self.re * rhs.re - self.im * rhs.im,
            im: self.re * rhs.im + self.im * rhs.re,
        }
    }

    pub fn scale(self, factor: i128) -> Self {
        Self {
            re: self.re * factor,
            im: self.im * factor,
        }
    }

    pub fn pow(self, mut exponent: u32) -> Self {
        let mut base = self;
        let mut acc = Self::one();
        while exponent > 0 {
            if exponent & 1 == 1 {
                acc = acc.mul(base);
            }
            exponent >>= 1;
            if exponent > 0 {
                base = base.mul(base);
            }
        }
        acc
    }

    pub fn div_exact(self, rhs: Self) -> Option<Self> {
        let norm = rhs.norm();
        if norm == 0 {
            return None;
        }
        let numerator = self.mul(rhs.conjugate());
        if numerator.re % norm == 0 && numerator.im % norm == 0 {
            Some(Self {
                re: numerator.re / norm,
                im: numerator.im / norm,
            })
        } else {
            None
        }
    }
}

pub fn units() -> [Gaussian; 4] {
    [
        Gaussian::new(1, 0),
        Gaussian::new(0, 1),
        Gaussian::new(-1, 0),
        Gaussian::new(0, -1),
    ]
}

#[cfg(test)]
mod tests {
    use super::Gaussian;

    #[test]
    fn gaussian_multiplication_and_norm() {
        let z = Gaussian::new(2, 1);
        let w = Gaussian::new(3, -4);
        assert_eq!(z.mul(w), Gaussian::new(10, -5));
        assert_eq!(z.norm(), 5);
        assert_eq!(w.norm(), 25);
        assert_eq!(z.mul(z.conjugate()), Gaussian::new(5, 0));
    }

    #[test]
    fn exact_gaussian_division() {
        let pi = Gaussian::new(2, 1);
        let z = pi.mul(Gaussian::new(7, -3));
        assert_eq!(z.div_exact(pi), Some(Gaussian::new(7, -3)));
        assert_eq!(Gaussian::new(1, 1).div_exact(pi), None);
    }
}
