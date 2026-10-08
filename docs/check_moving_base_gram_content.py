"""Checks for the moving-base Gram formula and exact weighted content."""

from fractions import Fraction as F
from functools import reduce
from itertools import combinations
from math import gcd, lcm


def gcd_all(values):
    return reduce(gcd, (abs(v) for v in values), 0)


def primitive(values):
    den = lcm(*(v.denominator for v in values))
    nums = [int(v * den) for v in values]
    common = gcd_all(nums)
    return [n // common for n in nums]


def weighted(nodes):
    rows = []
    for i, t in enumerate(nodes):
        d = F(1)
        for j, r in enumerate(nodes):
            if i != j:
                d *= t - r
        rows.append([t**k / d for k in range(5)])
    return rows


def check_pair_gcd():
    cases = 0
    for L in range(1, 31):
        for A in range(-80, 81):
            for B in range(A + 1, 81):
                delta = B - A
                if (A * A + L * L) % delta:
                    continue
                C = (A * A + L * L) // delta
                d = gcd_all((delta, A, C))
                assert L % d == 0, (A, B, L, d)
                G = (delta * delta, delta * A, A * A + L * L)
                assert gcd_all(G) == delta * d
                cases += 1
    return cases


def check_weighted_content(Xs, L, pair=(0, 1)):
    A, B = (Xs[i] for i in pair)
    delta = B - A
    assert delta > 0 and (A * A + L * L) % delta == 0
    C = (A * A + L * L) // delta
    d = gcd_all((delta, A, C))
    assert L % d == 0
    a, b, c = delta // d, A // d, C // d

    # Actual central vector with one point at infinity.
    ts = [F(x, L) for x in Xs]
    actual = []
    for i, t in enumerate(ts):
        D = reduce(lambda u, r: u * (t - r),
                   (r for j, r in enumerate(ts) if i != j), F(1))
        actual.append((t * t + 1) ** 2 / D)
    actual.append(F(-1))
    actual = primitive(actual)

    # Normalized finite nodes and Q0 evaluation.
    normalized = [F(x - A, delta) for x in Xs]
    finite = [normalized[i] for i in range(len(normalized))
              if i not in pair]
    # The selected rows are 0 and 1, followed by the three others.
    T = [F(0), F(1)] + finite
    assert len(T) == 5 and len(set(T)) == 5
    P = lambda t: a * t * t + 2 * b * t + c
    normalized_raw = []
    for i, t in enumerate(T):
        D = reduce(lambda u, r: u * (t - r),
                   (r for j, r in enumerate(T) if i != j), F(1))
        normalized_raw.append(P(t) ** 2 / D)
    normalized_raw.append(F(-a * a))
    assert actual == primitive(normalized_raw) or actual == [
        -v for v in primitive(normalized_raw)
    ]

    E = weighted(T)
    ell = lcm(*(entry.denominator for row in E for entry in row))
    N = [[int(ell * entry) for entry in row] for row in E]
    Nhat = N + [[0, 0, 0, 0, -ell]]
    f = [c * c, 4 * b * c, 2 * a * c + 4 * b * b, 4 * a * b, a * a]
    values = [sum(row[k] * f[k] for k in range(5)) for row in Nhat]
    content = gcd_all(values)
    assert primitive(normalized_raw) == [v // content for v in values]
    assert max(abs(v) for v in primitive(normalized_raw)) == (
        max(abs(v) for v in values) // content
    )


def main():
    pair_cases = check_pair_gcd()
    check_weighted_content((6, 10, 14, 22, 26), 4)
    check_weighted_content((2, 3, 5, 12, 17), 1)
    check_weighted_content((1, 2, 3, 4, 5), 1)
    print(f"PASS: {pair_cases} pair Gram reductions, 3 exact weighted-content checks")


if __name__ == "__main__":
    main()
