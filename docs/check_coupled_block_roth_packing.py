"""Exact finite checks for coupled_block_roth_packing.md.

Uses only the Python standard library. This audits Gaussian cancellation,
fixed power roots, and the primal/dual balanced-packing formulas. It does
not test Roth's theorem numerically. The general integral-vector
classification used by the LP dual is proved in balanced_pair_products.md.
"""

from fractions import Fraction
from itertools import combinations, product
from math import comb


ONE = (1, 0)
UNITS = ((1, 0), (0, 1), (-1, 0), (0, -1))
PRIMES = ((2, 1), (3, 2), (4, 1))
RATIONAL_PRIMES = (5, 13, 17)


def mul(z, w):
    return (z[0] * w[0] - z[1] * w[1], z[0] * w[1] + z[1] * w[0])


def conj(z):
    return (z[0], -z[1])


def norm(z):
    return z[0] ** 2 + z[1] ** 2


def power(z, exponent):
    assert exponent >= 0
    out = ONE
    while exponent:
        if exponent % 2:
            out = mul(out, z)
        z = mul(z, z)
        exponent //= 2
    return out


def nearest_integer(numerator, denominator):
    quotient, remainder = divmod(numerator, denominator)
    return quotient + (2 * remainder >= denominator)


def gaussian_gcd(z, w):
    while w != (0, 0):
        numerator = mul(z, conj(w))
        denominator = norm(w)
        quotient = tuple(nearest_integer(a, denominator) for a in numerator)
        multiple = mul(quotient, w)
        remainder = (z[0] - multiple[0], z[1] - multiple[1])
        assert norm(remainder) < norm(w)
        z, w = w, remainder
    return z


def signed_prime_product(exponents):
    out = ONE
    for prime, exponent in zip(PRIMES, exponents):
        out = mul(out, power(prime if exponent >= 0 else conj(prime), abs(exponent)))
    return out


def reduction(exponents, unit_indices, vector):
    """Return the exact primitive product and audit raw A = unit * c * g."""
    blocks = [
        mul(UNITS[unit], signed_prime_product(row))
        for unit, row in zip(unit_indices, exponents)
    ]
    assert all(norm(gaussian_gcd(block, conj(block))) == 1 for block in blocks)
    raw = ONE
    raw_unit = ONE
    for block, unit, coefficient in zip(blocks, unit_indices, vector):
        raw = mul(raw, power(block if coefficient >= 0 else conj(block), abs(coefficient)))
        raw_unit = mul(raw_unit, power(UNITS[unit if coefficient >= 0 else (-unit) % 4], abs(coefficient)))

    signed = [sum(v * row[p] for v, row in zip(vector, exponents)) for p in range(3)]
    positive = [sum(max(v * row[p], 0) for v, row in zip(vector, exponents)) for p in range(3)]
    negative = [sum(max(-v * row[p], 0) for v, row in zip(vector, exponents)) for p in range(3)]
    primitive = signed_prime_product(signed)
    content = 1
    raw_norm = 1
    primitive_norm = 1
    for p, e, plus, minus in zip(RATIONAL_PRIMES, signed, positive, negative):
        assert plus - minus == e
        assert plus + minus - abs(e) == 2 * min(plus, minus)
        content *= p ** min(plus, minus)
        raw_norm *= p ** (plus + minus)
        primitive_norm *= p ** abs(e)

    assert raw == mul(raw_unit, (content * primitive[0], content * primitive[1]))
    assert norm(raw) == raw_norm == content ** 2 * primitive_norm
    assert norm(primitive) == primitive_norm
    assert norm(gaussian_gcd(primitive, conj(primitive))) == 1
    assert norm(gaussian_gcd(raw, conj(raw))) == content ** 2
    assert (norm(primitive) == 1) == all(e == 0 for e in signed)
    return primitive, raw, raw_unit, content


def check_gaussian_cancellation():
    # Disjoint support; shared support in both orientations; repeated directions.
    families = (
        ((1, 0, 0), (0, 2, 0), (0, 0, -1)),
        ((1, 1, 0), (2, 0, -1), (0, -1, 1)),
        ((1, 0, 0), (1, 0, 0), (-1, 0, 0)),
    )
    units = ((0, 0, 0), (1, 2, 3), (3, 1, 2), (2, 3, 1))
    vectors = [v for v in product(range(-2, 3), repeat=3) if any(v)]
    cases = unit_cases = 0
    for exponents, unit_indices, vector in product(families, units, vectors):
        primitive, _, _, _ = reduction(exponents, unit_indices, vector)
        cases += 1
        unit_cases += norm(primitive) == 1

    # Opposite orientations yield a rational factor and a UNIT reduced product.
    primitive, raw, unit, content = reduction(((1, 0, 0), (-1, 0, 0)), (1, 0), (1, 1))
    assert (primitive, raw, unit, content) == (ONE, (0, 5), (0, 1), 5)
    # Equal directions in a quotient also collapse; v != 0 does not prevent this.
    primitive, raw, unit, content = reduction(((1, 0, 0), (1, 0, 0)), (0, 0), (1, -1))
    assert (primitive, raw, unit, content) == (ONE, (5, 0), ONE, 5)
    # A mixed example cancels 5^2 * 13^2 and retains bar(pi_13) * pi_17^4.
    primitive, _, _, content = reduction(families[1], (0, 0, 0), (2, -1, 3))
    assert content == 5 ** 2 * 13 ** 2
    assert primitive == mul(conj(PRIMES[1]), power(PRIMES[2], 4))
    return cases + 3, unit_cases + 2


def check_power_roots():
    base = ((1, 1, 0), (2, 0, -1), (0, -1, 1))
    vectors = ((1, -1, 0), (2, -1, 3), (0, 1, -2))
    cases = 0
    for ell, vector in product((2, 3, 4), vectors):
        scaled = tuple(tuple(ell * e for e in row) for row in base)
        gamma, _, _, _ = reduction(base, (1, 2, 3), vector)
        primitive, _, _, _ = reduction(scaled, (3, 1, 2), vector)
        assert norm(gamma) > 1
        assert primitive == power(gamma, ell)
        assert norm(primitive) == norm(gamma) ** ell
        cases += 1
    return cases


def check_balanced_packing():
    maxima = []
    for m in (4, 6, 8, 10):
        k = m // 2
        columns = tuple(combinations(range(m), k))
        pairs = tuple(combinations(range(m), 2))
        b = len(columns)
        a = Fraction(m - 2, 2 * (m - 1))
        weight = Fraction(1, k * (k - 1))
        pair_norms = [0] * len(pairs)
        for column in columns:
            signs = [1 if i in column else -1 for i in range(m)]
            values = [(signs[i] + signs[j]) // 2 for i, j in pairs]
            assert sum(weight * abs(v) for v in values) == 1
            for index, value in enumerate(values):
                pair_norms[index] += abs(value)
        assert all(Fraction(n, b) == a for n in pair_norms)

        # Uniform LP dual price 1/(a B). The prose classification reduces
        # every integral nonzero vector to spread D>=2 or these D=1 types.
        # Check the D>=2 lower bound already exceeds the dual threshold.
        assert Fraction(m * 2, 4 * (m - 1)) > a
        # For D=1, after a common translation n is binary with 2j ones.
        for j in range(1, k):
            total = 0
            probabilities = 0
            for x in range(2 * j + 1):
                if 0 <= k - x <= m - 2 * j:
                    count = comb(2 * j, x) * comb(m - 2 * j, k - x)
                    total += count * abs(x - j)
                    probabilities += count
            assert probabilities == b
            hypergeometric_norm = Fraction(total, b)
            closed_norm = Fraction(j * (k - j), k) * Fraction(comb(k, j) ** 2, comb(2 * k, 2 * j))
            assert hypergeometric_norm == closed_norm >= a
            assert (hypergeometric_norm == a) == (j in (1, k - 1))

        primal = len(pairs) * weight
        dual = b * Fraction(1, 1) / (a * b)
        assert primal == dual == Fraction(2 * (m - 1), m - 2) < 4
        maxima.append((m, primal))
    return maxima


def main():
    cases, unit_cases = check_gaussian_cancellation()
    roots = check_power_roots()
    maxima = check_balanced_packing()
    print("PASS: {} exact Gaussian cancellations ({} unit reductions); {} fixed power roots.".format(cases, unit_cases, roots))
    print("PASS: balanced packing primal/dual capacities: " + ", ".join("M={} -> {}".format(m, value) for m, value in maxima))
    print("Roth's theorem and the general integral-vector classification are prose inputs, not numerical tests.")


if __name__ == "__main__":
    main()
