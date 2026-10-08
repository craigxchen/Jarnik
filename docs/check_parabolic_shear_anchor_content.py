"""Exact checks for the parabolic shear's actual-anchor content charge."""

from fractions import Fraction
from itertools import product
from math import comb, gcd

from check_gaussian_reflection_replacement import gconj, gmul, gnorm
from check_mobius_conductor_transfer import least_norm
from check_mobius_source_aligned_compression import gaussian_valuation


def power(z, exponent):
    out = (1, 0)
    for _ in range(exponent):
        out = gmul(out, z)
    return out


def valuation(n, p):
    assert n
    n = abs(n)
    count = 0
    while n % p == 0:
        n //= p
        count += 1
    return count


def clamp(value, lower, upper):
    return max(lower, min(value, upper))


def width(profile, lower, upper):
    values = [clamp(v, lower, upper) for v in profile]
    return max(values) - min(values)


def local_audit():
    count = 0
    for e in range(1, 25):
        for anchor in range(e + 1):
            for t in range(1, 25):
                clipped = width([-anchor, 0, e - anchor], -t, t)
                assert clipped == min(t, anchor) + min(t, e - anchor)
                assert clipped <= t + min(anchor, e - anchor)
                if anchor in (0, e):
                    assert clipped <= t
                count += 1
    return count


def gaussian_audit():
    primes = [(2, 1), (3, 2), (4, 1), (5, 2)]
    count = sharp = 0
    # These are literal primitive physical circles. At the chosen anchor,
    # relative half-angle rows have signed depths e_i-e_j, with no factor 2.
    for pi, e, middle in product(primes, range(1, 5), range(5)):
        if middle > e:
            continue
        p = gnorm(pi)
        bar_pi = gconj(pi)
        physical = [power(bar_pi, e),
                    gmul(power(pi, middle), power(bar_pi, e - middle)),
                    power(pi, e)]
        if middle in (0, e):
            physical = [physical[0], physical[2]]
            anchor = 0 if middle == 0 else 1
        else:
            anchor = 1
        g_anchor = gcd(*physical[anchor])
        assert g_anchor == p ** min(middle, e - middle)
        levels = [0, e] if len(physical) == 2 else [0, middle, e]
        signed = [level - middle for level in levels]
        half_rows = [power(pi if s >= 0 else bar_pi, abs(s)) for s in signed]
        assert half_rows[anchor] == (1, 0)
        assert least_norm(half_rows) == p ** e
        for a in range(1, 6):
            for b in list(range(-8, 9)) + [p, p * p]:
                if b == 0 or gcd(a, b) != 1:
                    continue
                U, V = (2 * a, -b), (0, b)
                endpoints = [gaussian_valuation(V, pi) - gaussian_valuation(U, pi),
                             gaussian_valuation(gconj(U), pi)
                             - gaussian_valuation(gconj(V), pi)]
                lo, hi = sorted(endpoints)
                clipped = width(signed, lo, hi)
                bound = abs(b) * (b * b + 4 * a * a) * g_anchor
                assert valuation(bound, p) >= clipped
                if b % p == 0:
                    t = valuation(b, p)
                    assert (lo, hi) == (-t, t)
                elif (b * b + 4 * a * a) % p == 0:
                    assert hi - lo == valuation(b * b + 4 * a * a, p)
                else:
                    assert lo == hi == 0
                images = []
                for x, y in half_rows:
                    X, Y = a * x + b * y, a * y
                    content = gcd(X, Y)
                    assert a * a % content == 0
                    images.append((X // content, Y // content))
                N_prime = least_norm(images)
                assert N_prime % (p ** clipped) == 0
                count += 1
        if e % 2 == 0 and middle == e // 2:
            t = e // 2
            b = p ** t
            assert width(signed, -t, t) == 2 * t
            assert valuation(b, p) == t
            assert valuation(b * (b * b + 4) * g_anchor, p) == 2 * t
            sharp += 1
    return count, sharp


def constants(m):
    F = Fraction
    k = m // 4
    s = m - 2 * k
    A = F(k, s - 1)
    c = F(comb(k, 2)) + A * (F(k) - F(m, 2)) ** 2
    q = F(m) * (1 + 2 * A) / (8 * c)
    F_m = (m - 1) ** 2 // 4
    Q_m = m * m // 4
    ell = F(1, 4 * (m - 1)) + F(Q_m) * q / (m * (m - 1))
    g = F(m, 4 * F_m)
    h = ell + 2 * g
    def pair_exponent(n):
        return F(n * n // 4, n * (n - 1))
    u = min(pair_exponent(n) + pair_exponent(m - n + 1)
            for n in range(2, m))
    eta = 1 - q - F(3, 2) * u - 3 * h
    return q, u, h, eta


def exponent_audit():
    assert constants(128)[3] == Fraction(29828621, 231197785)
    assert constants(128)[3] > Fraction(1, 8)
    count = 0
    for m in (64, 96, 128, 192, 256, 512, 1024):
        q, u, h, eta = constants(m)
        assert eta > 0
        assert q <= Fraction(8, m)
        assert u <= Fraction(1, 2) + Fraction(2, m)
        assert h <= Fraction(7, m)
        assert eta >= Fraction(1, 4) - Fraction(32, m)
        if m >= 128:
            assert eta > Fraction(1, 8)
        assert abs(float(m * (Fraction(1, 4) - eta)) - 15.25) < 0.6
        count += 1
    return count


def condition_audit():
    count = 0
    # kappa is the larger root of x+1/x=2+b^2/a^2. Isolate it
    # between rational lower/upper certificates, without floating point.
    for a in range(1, 21):
        for b in range(1, 41):
            if gcd(a, b) != 1:
                continue
            T = Fraction(2 * a * a + b * b, a * a)
            # T-1 <= kappa <= T. The following inequalities hold already
            # at the lower bound, hence for the true condition number.
            lower = T - 1
            assert b * b <= a * a * lower
            assert b * b + 4 * a * a <= 5 * a * a * lower
            assert (b * (b * b + 4 * a * a)) ** 2 <= 25 * a ** 6 * lower ** 3
            count += 1
    return count


def main():
    local = local_audit()
    gaussian, sharp = gaussian_audit()
    exponents = exponent_audit()
    conditions = condition_audit()
    print(f"PASS: {local} exact clipped intervals; {gaussian} actual Gaussian "
          f"source/image checks; {sharp} necessary-anchor-content fixtures.")
    print(f"PASS: {exponents} exact exponent audits; {conditions} rational "
          "condition-number certificates.")


if __name__ == "__main__":
    main()
