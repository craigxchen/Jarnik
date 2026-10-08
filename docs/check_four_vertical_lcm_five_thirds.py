"""Exact checks for four_vertical_lcm_five_thirds.md; no floating point."""

from itertools import combinations, product
from fractions import Fraction
from math import prod

from check_endpoint_vertical_line_gcd import (
    exact_div,
    gcd_gaussian,
    gcd_many,
    lcm_many,
    norm,
)


def check_exponents() -> int:
    checks = 0
    for exponents in product(range(8), repeat=4):
        e1, e2, e3, e4 = sorted(exponents, reverse=True)
        pair_sum = sum(min(exponents[i], exponents[j])
                       for i, j in combinations(range(4), 2))
        assert 3 * e1 + pair_sum >= 2 * sum(exponents)
        assert 3 * e1 + pair_sum - 2 * sum(exponents) == e1 - e2 + e4
        checks += 1
    return checks


def check_vertical_tuples() -> int:
    checks = 0
    for d in range(1, 7):
        for heights in combinations(range(1, 17), 4):
            ws = [(-d, t) for t in heights]
            lnorm = norm(lcm_many(ws))
            norms_product = prod(map(norm, ws))
            pair_product = prod(norm(gcd_gaussian(ws[i], ws[j]))
                                for i, j in combinations(range(4), 2))
            triples = [gcd_many([ws[j] for j in subset])
                       for subset in combinations(range(4), 3)]
            gnorm = norm(gcd_many(ws))
            hnorm = norm(lcm_many(triples))
            snorm = prod(norm(exact_div(ws[j], gcd_gaussian(
                ws[j], lcm_many([ws[k] for k in range(4) if k != j]))))
                for j in range(4))

            assert lnorm**3 * pair_product >= norms_product**2
            assert lnorm**2 * hnorm * gnorm == norms_product * snorm
            assert lnorm * pair_product * gnorm == norms_product * prod(map(norm, triples))

            differences = prod(heights[j] - heights[i]
                               for i, j in combinations(range(4), 2))
            assert lnorm**3 * d**6 * differences >= norms_product**2
            t1, t2, t3, t4 = heights
            assert lnorm**3 * d**6 >= t1**4 * t2**3 * t3**2 * t4
            assert lnorm**3 * d**6 >= t1**10
            checks += 1
    return checks


def check_relaxed_obstruction() -> None:
    blocks = [frozenset(s) for size in (2, 3)
              for s in combinations(range(4), size)]
    factors = [{b for b in blocks if j in b} for j in range(4)]
    assert len(blocks) == 10
    assert all(len(f) == 6 for f in factors)
    assert all(len(factors[i] & factors[j]) == 3
               for i, j in combinations(range(4), 2))
    assert all(len(set.intersection(*(factors[j] for j in subset))) == 1
               for subset in combinations(range(4), 3))
    assert not set.intersection(*factors)
    assert len(set.union(*factors)) == 10


def check_general_cardinality() -> int:
    checks = 0
    for k in range(2, 31):
        r = k // 2
        denominator = r * (r + 1)
        for s in range(k + 1):
            assert denominator - 2 * r * s + s * (s - 1) == (s - r) * (s - r - 1)
            assert (s - r) * (s - r - 1) >= 0
            checks += 1
        height_numerators = [2 * r - j + 1 for j in range(1, k + 1)]
        assert min(height_numerators) >= 0
        assert sum(height_numerators) == r * (2 * r + 1)
        assert sum(height_numerators) >= k * (k - 1) // 2
    return checks


def check_varying_contents() -> int:
    checks = 0
    candidates = [(-d, t) for d in range(1, 4) for t in range(1, 6)]
    for unordered in combinations(candidates, 4):
        ws = sorted(unordered, key=lambda w: Fraction(w[1], -w[0]))
        ds = [-w[0] for w in ws]
        xs = [Fraction(w[1], -w[0]) for w in ws]
        if len(set(xs)) != 4:
            continue
        for i, j in combinations(range(4), 2):
            determinant = abs(ws[i][0] * ws[j][1] - ws[i][1] * ws[j][0])
            assert norm(gcd_gaussian(ws[i], ws[j])) <= determinant
            assert determinant == ds[i] * ds[j] * (xs[j] - xs[i])
        lnorm = norm(lcm_many(ws))
        # Sixth power of (11): |L|^6=N(L)^3.
        assert lnorm**3 >= prod(ds) * xs[0]**4 * xs[1]**3 * xs[2]**2 * xs[3]
        assert lnorm**3 >= xs[0]**10
        checks += 1
    return checks


if __name__ == "__main__":
    print(f"PASS: {check_exponents()} exact valuation checks")
    print(f"PASS: {check_vertical_tuples()} four-cofactor checks")
    check_relaxed_obstruction()
    print("PASS: ten-block relaxed obstruction")
    print(f"PASS: {check_general_cardinality()} general-cardinality multiplicity checks")
    print(f"PASS: {check_varying_contents()} varying-content four-cofactor checks")
