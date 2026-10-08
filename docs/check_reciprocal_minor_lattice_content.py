"""Exact audit for reciprocal_minor_lattice_content.md."""

from __future__ import annotations

from fractions import Fraction
from itertools import combinations
from math import gcd, lcm
from random import Random


Gaussian = tuple[int, int]


def gmul(a: Gaussian, b: Gaussian) -> Gaussian:
    return (a[0] * b[0] - a[1] * b[1],
            a[0] * b[1] + a[1] * b[0])


def gnorm(a: Gaussian) -> int:
    return a[0] * a[0] + a[1] * a[1]


def grem(a: Gaussian, b: Gaussian) -> Gaussian:
    """Remainder after nearest Gaussian-integer division."""
    den = gnorm(b)
    real_num = a[0] * b[0] + a[1] * b[1]
    imag_num = a[1] * b[0] - a[0] * b[1]
    qr = (2 * real_num + den) // (2 * den)
    qi = (2 * imag_num + den) // (2 * den)
    return (a[0] - qr * b[0] + qi * b[1],
            a[1] - qr * b[1] - qi * b[0])


def ggcd(a: Gaussian, b: Gaussian) -> Gaussian:
    while b != (0, 0):
        a, b = b, grem(a, b)
    return a


def gquot(a: Gaussian, b: Gaussian) -> Gaussian:
    den = gnorm(b)
    real_num = a[0] * b[0] + a[1] * b[1]
    imag_num = a[1] * b[0] - a[0] * b[1]
    assert real_num % den == 0 and imag_num % den == 0
    return (real_num // den, imag_num // den)


def glcm(a: Gaussian, b: Gaussian) -> Gaussian:
    return gquot(gmul(a, b), ggcd(a, b))


def det(a: Gaussian, b: Gaussian) -> int:
    return a[0] * b[1] - a[1] * b[0]


def differ_by_common_unit(xs: list[Gaussian], ys: list[Gaussian]) -> bool:
    units = [(1, 0), (-1, 0), (0, 1), (0, -1)]
    return any([gmul(unit, x) for x in xs] == ys for unit in units)


def data(d: int, ts: list[int]) -> tuple[int, list[Gaussian], list[int], int, int]:
    A = (1, 0)
    for t in ts:
        A = glcm(A, (-d, t))
    M = gnorm(A)
    qs = [gquot(A, (-d, t)) for t in ts]
    deltas = [abs(det(qs[i], qs[j]))
              for i, j in combinations(range(len(ts)), 2)]
    h = gcd(*deltas)
    denominators = []
    for i, j in combinations(range(len(ts)), 2):
        ni = d * d + ts[i] * ts[i]
        nj = d * d + ts[j] * ts[j]
        num = d * (ts[j] - ts[i])
        denominators.append(ni * nj // gcd(ni * nj, num))
    B = lcm(*denominators)
    return M, qs, deltas, h, B


def gcd_all(values: list[int]) -> int:
    result = 0
    for value in values:
        result = gcd(result, abs(value))
    return result


def audit_instance(d: int, ts: list[int]) -> None:
    M, qs, deltas, h, B = data(d, ts)

    # Exact lcm/content formula, with no odd-prime restriction.
    assert B == M // gcd(M, h)

    # Gaussian gcd one implies ordinary coordinate content one.
    common_gcd = (0, 0)
    for q in qs:
        common_gcd = ggcd(common_gcd, q)
    assert gnorm(common_gcd) == 1
    assert gcd_all([coordinate for q in qs for coordinate in q]) == 1

    # Maximal minors of [q_0 ... q_(k-1) M e_1 M e_2].
    augmented_minors = deltas[:]
    for x, y in qs:
        augmented_minors.extend([M * x, M * y])
    augmented_minors.append(M * M)
    assert gcd_all(augmented_minors) == gcd(M, h)

    # Check the determinant formula independently of Gaussian lcm arithmetic.
    expected = []
    for i, j in combinations(range(len(ts)), 2):
        ni = d * d + ts[i] * ts[i]
        nj = d * d + ts[j] * ts[j]
        numerator = M * d * (ts[j] - ts[i])
        assert numerator % (ni * nj) == 0
        expected.append(numerator // (ni * nj))
    assert deltas == expected


def audit_scaling(d: int, ts: list[int]) -> None:
    M, qs, _, h, B = data(d, ts)
    for c in (2, 3, 5, 12):
        Mc, qsc, _, hc, Bc = data(c * d, [c * t for t in ts])
        assert Mc == c * c * M
        # The Gaussian lcm may change by a unit, so all q's share one unit.
        assert differ_by_common_unit(qs, qsc)
        assert hc == h
        assert Bc == c * c * M // gcd(c * c * M, h)
        target = Fraction(ts[0] ** 4, d ** 2)
        target_c = Fraction((c * ts[0]) ** 4, (c * d) ** 2)
        assert target_c == c * c * target
        assert Bc / target_c <= B / target


def audit_family(limit: int = 151) -> int:
    checked = 0
    for d in range(1, limit + 1, 2):
        if gcd(d, 6) != 1:
            continue
        ts = [d + 4 * j for j in range(4)]
        assert gcd_all([d, *ts]) == 1
        us = [(2 * j, d + 2 * j) for j in range(4)]
        Ns = [gnorm(u) for u in us]
        for i, j in combinations(range(4), 2):
            assert gnorm(ggcd(us[i], us[j])) == 1

        explicit_A = (1, 1)
        for u in us:
            explicit_A = gmul(explicit_A, u)
        explicit_qs = []
        for j in range(4):
            q = (1, 0)
            for ell in range(4):
                if ell != j:
                    q = gmul(q, us[ell])
            explicit_qs.append(q)

        M, qs, deltas, h, B = data(d, ts)
        assert M == gnorm(explicit_A) == 2 * Ns[0] * Ns[1] * Ns[2] * Ns[3]
        assert differ_by_common_unit(explicit_qs, qs)
        expected_deltas = [
            2 * d * (j - i)
            * Ns[[ell for ell in range(4) if ell not in (i, j)][0]]
            * Ns[[ell for ell in range(4) if ell not in (i, j)][1]]
            for i, j in combinations(range(4), 2)
        ]
        assert deltas == expected_deltas
        assert h == 2 * d
        assert gcd(M, h) == 2 * d
        assert B == M // (2 * d) == Ns[0] * Ns[1] * Ns[2] * Ns[3] // d

        common_gcd = (0, 0)
        for q in explicit_qs:
            common_gcd = ggcd(common_gcd, q)
        assert gnorm(common_gcd) == 1
        checked += 1
    return checked


def main() -> None:
    fixed = [
        (2, [6, 8, 16, 42]),
        (1, [3, 5, 13, 31]),
        (35, [931, 2376]),       # two-point content has no 2d bound
        (75, [221, 261, 289, 325]),
    ]
    for d, ts in fixed:
        audit_instance(d, ts)
        audit_scaling(d, ts)

    rng = Random(20260913)
    random_trials = 0
    for k in (2, 3, 4, 5, 6):
        for _ in range(200):
            d = rng.randrange(1, 80)
            ts = sorted(rng.sample(range(1, 500), k))
            audit_instance(d, ts)
            random_trials += 1

    family_cases = audit_family()
    print(
        "PASS: exact content identity on "
        f"{len(fixed) + random_trials} instances; "
        f"scaling on {len(fixed)} instances; family on {family_cases} values of d"
    )


if __name__ == "__main__":
    main()
