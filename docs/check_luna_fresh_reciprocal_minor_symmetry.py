"""Exact finite audit for luna_fresh_reciprocal_minor_symmetry.md."""

from __future__ import annotations

from math import gcd, lcm
from random import Random


Gaussian = tuple[int, int]


def gdiv(a: Gaussian, b: Gaussian) -> Gaussian:
    """Remainder after nearest Gaussian division."""
    den = b[0] * b[0] + b[1] * b[1]
    x_num = a[0] * b[0] + a[1] * b[1]
    y_num = a[1] * b[0] - a[0] * b[1]
    q = ((2 * x_num + den) // (2 * den),
         (2 * y_num + den) // (2 * den))
    return (
        a[0] - q[0] * b[0] + q[1] * b[1],
        a[1] - q[0] * b[1] - q[1] * b[0],
    )


def gquot(a: Gaussian, b: Gaussian) -> Gaussian:
    """Exact quotient when ``b`` divides ``a`` in ``Z[i]``."""
    den = b[0] * b[0] + b[1] * b[1]
    real = (a[0] * b[0] + a[1] * b[1])
    imag = (a[1] * b[0] - a[0] * b[1])
    assert real % den == 0 and imag % den == 0, (a, b)
    return (real // den, imag // den)


def ggcd(a: Gaussian, b: Gaussian) -> Gaussian:
    while b != (0, 0):
        a, b = b, gdiv(a, b)
    return a


def glcm(a: Gaussian, b: Gaussian) -> Gaussian:
    g = ggcd(a, b)
    product = (a[0] * b[0] - a[1] * b[1], a[0] * b[1] + a[1] * b[0])
    return gquot(product, g)


def minor_lcm(d: int, ts: list[int]) -> int:
    value = 1
    for i in range(len(ts)):
        ni = d * d + ts[i] * ts[i]
        for j in range(i + 1, len(ts)):
            nj = d * d + ts[j] * ts[j]
            denominator = ni * nj
            numerator = d * (ts[j] - ts[i])
            value = lcm(value, denominator // gcd(denominator, numerator))
    return value


def lcm_norm(d: int, ts: list[int]) -> int:
    a: Gaussian = (1, 0)
    for t in ts:
        a = glcm(a, (-d, t))
    return a[0] * a[0] + a[1] * a[1]


def main() -> None:
    exact_examples = [
        (2, [6, 8, 16, 42], 8840),
        (1, [3, 5, 13, 31], 81770),
        (2, [6, 8, 10, 36], 44200),
    ]
    for d, ts, expected in exact_examples:
        b = minor_lcm(d, ts)
        n_a = lcm_norm(d, ts)
        assert b == expected and n_a == expected
        assert n_a % b == 0

    rng = Random(20260913)
    trials = 0
    for d in range(1, 31):
        for _ in range(1000):
            ts = sorted(rng.sample(range(d, d + 500), 4))
            b = minor_lcm(d, ts)
            n_a = lcm_norm(d, ts)
            assert n_a % b == 0, (d, ts, b, n_a)
            trials += 1
    print(f"verified B_4 | N(A) on {trials} exact random quadruples")


if __name__ == "__main__":
    main()
