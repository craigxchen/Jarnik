"""Exact checks of the source-norm sparsity reduction for endpoint maps."""

from fractions import Fraction
from math import gcd, prod

from check_endpoint_map_small_parameter_reconstruction import constants, reanchor
from check_gaussian_reflection_replacement import ggcd, gnorm
from check_mobius_conductor_transfer import fixtures
from check_mobius_reciprocal_stretch_grid import clips, gram, realize


def actual_divisibility():
    count = 0
    for rows in fixtures() + [[(1, 0), (1, 1), (2, 1), (3, 1)]]:
        source, N = realize(rows)
        for j, anchor in enumerate(source):
            anchored = reanchor(rows, j)
            for a in range(1, 4):
                for d in range(1, 4):
                    for b in range(-5, 6):
                        if gcd(a, b, d) != 1 or (a == d and b == 0):
                            continue
                        T, W = gram(a, b, d)
                        K = (b * b + (a + d) ** 2) * (b * b + (a - d) ** 2)
                        assert K == gnorm(W) == T * T - 4 * a * a * d * d > 0
                        Q = prod(p ** e for p, e in clips(a, b, d, anchor, N).items())
                        E = N // Q
                        assert N % Q == K % Q == (E * K) % N == 0
                        target, Nprime = realize([(a * x + b * y, d * y)
                                                  for x, y in anchored])
                        inverse_Q = prod(p ** e for p, e in
                                         clips(d, -b, a, target[j], Nprime).items())
                        assert inverse_Q == Q
                        assert ((Nprime // Q) * K) % Nprime == 0
                        count += 1
    return count


def bounds():
    for m in range(128, 10001):
        q, h, u, rho, log2D = constants(m)
        assert q <= Fraction(8, m)
        assert h <= Fraction(7, m)
        assert u <= Fraction(1, 2) + Fraction(2, m)
        assert log2D <= 14
        Aexponent = h / rho
        Bexponent = (u / 2 + h) / rho
        assert Aexponent <= Fraction(7, m - 4)
        assert Bexponent <= Fraction(1, 4) + Fraction(9, m - 4)
        assert log2D + 4 * Aexponent <= 15
        # log2(sqrt(5)) < 5/4, an exact certificate since 5^2<2^5.
        assert Fraction(5, 4) + u / 2 + log2D + 4 * Bexponent <= 17
        beta = Fraction(1, 4) + Fraction(8, m) + Fraction(23, m - 4)
        assert q + 2 * Aexponent + Bexponent <= beta
        assert 1 + Fraction(8, m) + Fraction(36, m - 4) < 2
    assert Fraction(1, 4) + Fraction(8, 128) + Fraction(23, 124) == Fraction(247, 496)
    assert Fraction(247, 496) < Fraction(1, 2)
    q, h, u, rho, _ = constants(68)
    exact68 = q + (u / 2 + 3 * h) / rho
    assert exact68 == Fraction(96152107, 195629280) < Fraction(1, 2)
    q, h, u, rho, _ = constants(100000)
    exact = q + (u / 2 + 3 * h) / rho
    assert abs(100000 * (exact - Fraction(1, 4)) - Fraction(61, 4)) < Fraction(1, 100)
    return 10001 - 128


def actual_norm_family():
    for n in range(1, 501):
        lower, upper = (2 * n, -1), (2 * n, 1)
        N = 4 * n * n + 1
        assert gnorm(lower) == gnorm(upper) == N
        assert gnorm(ggcd(lower, upper)) == 1
        assert N <= 5 * n ** 4 < 16 * n ** 4
    return 500


def main():
    cases = actual_divisibility()
    exponents = bounds()
    norms = actual_norm_family()
    print(f"PASS: {cases} actual forward/inverse source-divisor checks.")
    print(f"PASS: {exponents} uniform exponent/constant checks; beta_128=247/496;")
    print("      exact beta_68<1/2 and asymptotic correction 61/(4m).")
    print(f"PASS: {norms} primitive actual endpoint-source norms 4n^2+1.")


if __name__ == '__main__':
    main()
