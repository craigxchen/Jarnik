#!/usr/bin/env python3
"""Exact finite audits for monomial_vector_subspace_height_barrier.md.

No approximation theorem or existence of endpoint families is tested.
"""

from fractions import Fraction
from itertools import combinations, product
from random import Random

from check_coupled_block_roth_packing import (
    ONE, UNITS, PRIMES, RATIONAL_PRIMES,
    mul, conj, norm, power, gaussian_gcd, signed_prime_product,
)


def nonnegative_prime_product(plus, minus):
    answer = ONE
    for prime, a, b in zip(PRIMES, plus, minus):
        assert a >= 0 and b >= 0
        answer = mul(answer, power(prime, a))
        answer = mul(answer, power(conj(prime), b))
    return answer


def check_height_case(exponents, units, vectors):
    assert not any(vectors[0])
    blocks = [mul(UNITS[u], signed_prime_product(row))
              for row, u in zip(exponents, units)]
    signed = [tuple(sum(v * row[p] for v, row in zip(vector, exponents))
                    for p in range(len(PRIMES))) for vector in vectors]
    lower = tuple(min(row[p] for row in signed) for p in range(len(PRIMES)))
    upper = tuple(max(row[p] for row in signed) for p in range(len(PRIMES)))
    denominator = nonnegative_prime_product(tuple(-a for a in lower), upper)
    common_norm = 1
    for p, lo, hi in zip(RATIONAL_PRIMES, lower, upper):
        common_norm *= p ** (hi - lo)
    assert norm(denominator) == common_norm
    coordinates = []
    for vector, row in zip(vectors, signed):
        raw = ONE
        raw_unit = ONE
        for v, block, u in zip(vector, blocks, units):
            raw = mul(raw, power(block if v >= 0 else conj(block), abs(v)))
            raw_unit = mul(raw_unit, power(UNITS[u if v >= 0 else -u % 4], abs(v)))
        coordinate = nonnegative_prime_product(
            tuple(e - lo for e, lo in zip(row, lower)),
            tuple(hi - e for e, hi in zip(row, upper)),
        )
        # chi_v carries the square of the raw Gaussian unit.
        coordinate = mul(power(raw_unit, 2), coordinate)
        assert mul(coordinate, conj(raw)) == mul(denominator, raw)
        assert norm(coordinate) == common_norm
        coordinates.append(coordinate)
    gcd_value = coordinates[0]
    for coordinate in coordinates[1:]:
        gcd_value = gaussian_gcd(gcd_value, coordinate)
    assert norm(gcd_value) == 1
    return common_norm == 1


def check_exact_heights():
    rng = Random(20260915)
    total = unit_vectors = 0
    families = (
        ((1, 0, 0), (0, 2, 0), (0, 0, -1)),
        ((1, 1, 0), (2, 0, -1), (0, -1, 1)),
        ((1, 0, 0), (1, 0, 0), (-1, 0, 0)),
        ((0, 0, 0), (0, 0, 0), (0, 0, 0)),
    )
    for family in families:
        for _ in range(50):
            q = rng.randrange(1, 6)
            vectors = [(0, 0, 0)] + [tuple(rng.randrange(-3, 4) for _ in range(3))
                                     for _ in range(q)]
            units = tuple(rng.randrange(4) for _ in range(3))
            unit_vectors += check_height_case(family, units, vectors)
            total += 1
    # A full collapse despite distinct exponent characters.
    unit_vectors += check_height_case(families[2], (1, 2, 3),
                                     [(0, 0, 0), (1, -1, 0), (1, 0, 1)])
    total += 1
    return total, unit_vectors


def check_range_identities():
    total = 0
    for a, b in product(range(-20, 21), repeat=2):
        assert 2 * (max(0, a, b) - min(0, a, b)) == abs(a) + abs(b) + abs(a-b)
        for shift in (-17, 0, 23):
            assert max(shift, shift+a, shift+b) - min(shift, shift+a, shift+b) == max(0, a, b) - min(0, a, b)
        total += 1
    return total


def balanced_vector(profile, cuts):
    assert sum(profile) % 2 == 0
    half = sum(profile) // 2
    return tuple(sum(profile[i] for i in cut) - half for cut in cuts)


def check_balanced_ranges():
    rng = Random(19062026)
    rows = []
    cases = 0
    for m in (4, 6, 8, 10, 12):
        cuts = tuple(combinations(range(m), m // 2))
        minimum = Fraction(m - 2, 2 * (m - 1))
        assert minimum > Fraction(1, 4)
        assert Fraction(3, 2) * minimum >= Fraction(1, 2)
        pool = {(0,) * len(cuts)}
        for _ in range(70):
            profile = [rng.randrange(-3, 4) for _ in range(m)]
            if sum(profile) % 2:
                profile[-1] += 1
            value = balanced_vector(profile, cuts)
            pool.add(value)
            if any(value):
                assert Fraction(sum(map(abs, value)), len(cuts)) >= minimum
        pool = sorted(pool)
        for q in range(1, 10):
            selected = rng.sample(pool, q + 1)
            width = Fraction(sum(max(col) - min(col) for col in zip(*selected)), len(cuts))
            bound = minimum if q == 1 else Fraction(3, 2) * minimum
            assert width >= bound > Fraction(q, 2 * (q+1))
            # A common arbitrary integral shift never alters projective width.
            shift = rng.choice(pool)
            translated = [tuple(a+b for a, b in zip(value, shift)) for value in selected]
            assert width == Fraction(sum(max(col) - min(col) for col in zip(*translated)), len(cuts))
            cases += 1
        rows.append((m, minimum, Fraction(3, 2) * minimum))
    # Exact three-character attaining configuration at M=6.
    cuts = tuple(combinations(range(6), 3))
    a = balanced_vector((1, 1, 0, 0, 0, 0), cuts)
    b = tuple(-x for x in balanced_vector((0, 0, 1, 1, 0, 0), cuts))
    zero = (0,) * len(cuts)
    width = Fraction(sum(max(col) - min(col) for col in zip(zero, a, b)), len(cuts))
    assert width == Fraction(3, 5)
    return cases + 1, rows


if __name__ == "__main__":
    total, unit_vectors = check_exact_heights()
    ranges = check_range_identities()
    balanced, bounds = check_balanced_ranges()
    print(f"PASS: {total} exact primitive projective height cases ({unit_vectors} unit-height vectors).")
    print(f"PASS: {ranges} range identities with translations; {balanced} balanced exponent-set checks.")
    print("Scalar and three-point range bounds: " + ", ".join(f"M={m}: {a}, {b}" for m, a, b in bounds))
    print("The Subspace theorem and all-exponent balanced minimum are prose inputs; no endpoint realization is asserted.")
