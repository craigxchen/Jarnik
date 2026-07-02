use std::cmp::Ordering;

pub fn isqrt_u128(n: u128) -> u128 {
    if n < 2 {
        return n;
    }

    let mut lo = 1_u128;
    let mut hi = 1_u128 << ((128 - n.leading_zeros() as usize + 1) / 2);
    while lo <= hi {
        let mid = lo + (hi - lo) / 2;
        match mid.checked_mul(mid).map(|sq| sq.cmp(&n)) {
            Some(Ordering::Equal) => return mid,
            Some(Ordering::Less) => lo = mid + 1,
            _ => hi = mid - 1,
        }
    }
    hi
}

pub fn square_root_i128(n: i128) -> Option<i128> {
    if n < 0 {
        return None;
    }
    let root = isqrt_u128(n as u128);
    if root * root == n as u128 {
        Some(root as i128)
    } else {
        None
    }
}

pub fn is_square_i128(n: i128) -> bool {
    square_root_i128(n).is_some()
}

pub fn mod_i128(n: i128, modulus: u128) -> u128 {
    let m = modulus as i128;
    let mut r = n % m;
    if r < 0 {
        r += m;
    }
    r as u128
}

pub fn gcd_i128(mut a: i128, mut b: i128) -> i128 {
    a = a.abs();
    b = b.abs();
    while b != 0 {
        let r = a % b;
        a = b;
        b = r;
    }
    a
}

pub fn mod_inverse_i128(a: i128, modulus: i128) -> Option<i128> {
    if modulus == 1 {
        return Some(0);
    }

    let (mut old_r, mut r) = (a, modulus);
    let (mut old_s, mut s) = (1_i128, 0_i128);

    while r != 0 {
        let quotient = old_r / r;
        let next_r = old_r - quotient * r;
        old_r = r;
        r = next_r;

        let next_s = old_s - quotient * s;
        old_s = s;
        s = next_s;
    }

    if old_r.abs() != 1 {
        return None;
    }
    let mut inverse = old_s % modulus;
    if inverse < 0 {
        inverse += modulus;
    }
    Some(inverse)
}

pub fn div_exact(n: i128, d: i128) -> Option<i128> {
    if d != 0 && n % d == 0 {
        Some(n / d)
    } else {
        None
    }
}

pub fn parse_i128(s: &str) -> Result<i128, String> {
    s.parse::<i128>()
        .map_err(|err| format!("could not parse `{}` as an integer: {}", s, err))
}

pub fn parse_u64_list(s: &str) -> Result<Vec<u64>, String> {
    if s.trim().is_empty() {
        return Ok(Vec::new());
    }
    s.split(',')
        .map(|part| {
            part.trim()
                .parse::<u64>()
                .map_err(|err| format!("could not parse prime `{}`: {}", part, err))
        })
        .collect()
}

pub fn cmp_angle(a: &(i128, i128), b: &(i128, i128)) -> Ordering {
    fn upper_half((x, y): &(i128, i128)) -> bool {
        *y > 0 || (*y == 0 && *x >= 0)
    }

    let a_upper = upper_half(a);
    let b_upper = upper_half(b);
    if a_upper != b_upper {
        return b_upper.cmp(&a_upper);
    }

    let cross = a.0 * b.1 - a.1 * b.0;
    if cross > 0 {
        Ordering::Less
    } else if cross < 0 {
        Ordering::Greater
    } else {
        a.cmp(b)
    }
}

#[cfg(test)]
mod tests {
    use super::{isqrt_u128, square_root_i128};

    #[test]
    fn exact_square_roots() {
        assert_eq!(isqrt_u128(0), 0);
        assert_eq!(isqrt_u128(1), 1);
        assert_eq!(isqrt_u128(2), 1);
        assert_eq!(isqrt_u128(15), 3);
        assert_eq!(isqrt_u128(16), 4);
        assert_eq!(square_root_i128(144), Some(12));
        assert_eq!(square_root_i128(145), None);
        assert_eq!(square_root_i128(-1), None);
    }
}
