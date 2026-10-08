"""Exact checks for the moving-base good-cut content table."""

from fractions import Fraction as F
from functools import reduce
from itertools import combinations
from math import gcd, lcm


P = 13


def valuation(n, p=P):
    """Return v_p(n) for a nonzero integer."""
    assert n
    n = abs(n)
    value = 0
    while n % p == 0:
        n //= p
        value += 1
    return value


def gcd_all(values):
    return reduce(gcd, (abs(value) for value in values), 0)


def hensel_root_minus_one(e):
    """The lift modulo 13^e of the root 5 modulo 13."""
    assert e >= 1
    root = 5
    modulus = P
    for _ in range(1, e):
        new_modulus = modulus * P
        lifts = [
            root + digit * modulus
            for digit in range(P)
            if ((root + digit * modulus) ** 2 + 1) % new_modulus == 0
        ]
        assert len(lifts) == 1
        root = lifts[0]
        modulus = new_modulus
    assert 0 <= root < modulus
    assert (root * root + 1) % modulus == 0
    return root


def clustered_cotangents(e):
    modulus = P**e
    rho = hensel_root_minus_one(e)

    # Exactly one shift modulo 13 lifts the congruence by another digit.
    bad = [
        shift
        for shift in range(P)
        if valuation((rho + shift * modulus) ** 2 + 1) > e
    ]
    assert len(bad) == 1
    bad = bad[0]

    start = next(
        shift
        for shift in range(P)
        if all((shift + offset) % P != bad for offset in range(3))
    )
    r = rho + start * modulus
    cluster = [r + offset * modulus for offset in range(3)]
    assert all(valuation(value * value + 1) == e for value in cluster)
    return [1, 2, *cluster]


def weighted_rows(nodes):
    rows = []
    for i, node in enumerate(nodes):
        denominator = F(1)
        for j, other in enumerate(nodes):
            if i != j:
                denominator *= node - other
        rows.append([node**degree / denominator for degree in range(5)])
    return rows


def check_abstract_table():
    expected = {
        (0, 1): 0,
        (0, 2): 1,
        (0, 3): 2,
        (0, 4): 2,
        (1, 2): 2,
        (1, 3): 2,
        (1, 4): 1,
        (1, 5): 0,
    }
    ell_total = 0
    content_total = 0

    labels = range(1, 6)
    for size in range(1, 6):
        for subset_tuple in combinations(labels, size):
            subset = set(subset_tuple)
            h = int({1, 2} <= subset)
            ell_order = 0 if h else size - 1
            raw_orders = [2 * h]
            raw_orders.extend(
                2 * h + (3 - size) * int(label in subset)
                for label in labels
            )
            content_order = ell_order + min(raw_orders)
            assert content_order == expected[h, size]
            ell_total += ell_order
            content_total += content_order

    assert ell_total == 29
    assert content_total == 38


def check_family(e):
    slopes = clustered_cotangents(e)
    edge_cotangents = [
        F(slopes[i] * slopes[j] + 1, slopes[j] - slopes[i])
        for i, j in combinations(range(5), 2)
    ]
    common_denominator = lcm(*(value.denominator for value in edge_cotangents))
    assert valuation(common_denominator) == 0

    xs = [common_denominator * slope for slope in slopes]
    a_value, b_value = xs[:2]
    delta = b_value - a_value
    quotient = (a_value * a_value + common_denominator**2) // delta
    d = gcd_all((delta, a_value, quotient))
    gram = (delta // d, a_value // d, quotient // d)
    assert gram == (1, 1, 2)
    assert gram[0] * gram[2] - gram[1] ** 2 == 1

    nodes = [F(value - a_value, delta) for value in xs]
    assert nodes == [F(0), F(1), *map(lambda value: F(value - 1), slopes[2:])]

    rows = weighted_rows(nodes)
    ell = lcm(*(entry.denominator for row in rows for entry in row))
    assert valuation(ell) == 2 * e

    aa, bb, cc = gram
    coefficients = [
        cc * cc,
        4 * bb * cc,
        2 * aa * cc + 4 * bb * bb,
        4 * aa * bb,
        aa * aa,
    ]
    integer_rows = [[int(ell * entry) for entry in row] for row in rows]
    integer_rows.append([0, 0, 0, 0, -ell])
    evaluations = [
        sum(row[index] * coefficients[index] for index in range(5))
        for row in integer_rows
    ]
    content = gcd_all(evaluations)
    assert valuation(content) == 2 * e

    # The cleared rational pair cotangents are integral in the X/L chart.
    for i, j in combinations(range(5), 2):
        numerator = xs[i] * xs[j] + common_denominator**2
        denominator = xs[j] - xs[i]
        assert numerator % denominator == 0

    return common_denominator, max(abs(value) for value in evaluations) // content


def main():
    check_abstract_table()
    results = [check_family(e) for e in range(1, 9)]
    print("PASS: good-cut table totals ell=29 and C_eval=38")
    print("PASS: fixed Q0 determinant-one family for e=1,...,8")
    print("least pair-cotangent denominators and primitive heights:")
    for e, (denominator, height) in enumerate(results, 1):
        print(f"  e={e}: L={denominator}, U={height}")


if __name__ == "__main__":
    main()
