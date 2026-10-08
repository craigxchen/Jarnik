"""Exact Fourier ideal-lattice identities for literal Gaussian blocks.

No fixture is asserted to have endpoint-small source arguments. Minimum
classification is proved by the explicit orthogonal operator, not search.
"""

from fractions import Fraction
from itertools import product
from math import prod
from random import Random

from check_five_row_gradient_cofactor_stress import conj, mul, norm, power, split_primes


ZERO, ONE = (0, 0), (1, 0)
UNITS = (ONE, (0, 1), (-1, 0), (0, -1))


def add(a, b):
    return a[0] + b[0], a[1] + b[1]


def sub(a, b):
    return a[0] - b[0], a[1] - b[1]


def scale(a, n):
    return a[0] * n, a[1] * n


def gsum(values):
    result = ZERO
    for value in values:
        result = add(result, value)
    return result


def gdiv(a, b):
    numerator = mul(a, conj(b))
    denominator = norm(b)
    return Fraction(numerator[0], denominator), Fraction(numerator[1], denominator)


def divides(divisor, value):
    numerator = mul(value, conj(divisor))
    denominator = norm(divisor)
    return numerator[0] % denominator == numerator[1] % denominator == 0


def character(v, x):
    return -1 if bin(v & x).count("1") % 2 else 1


def delta(order, label=0, value=ONE):
    result = [ZERO] * order
    result[label] = value
    return result


def translate(vector, v):
    return [vector[y ^ v] for y in range(len(vector))]


def apply_block(vector, v, block, adjoint=False):
    a, b = block
    if adjoint:
        b = -b
    return [add(scale(vector[y], a), mul((0, b), vector[y ^ v]))
            for y in range(len(vector))]


def apply_all(vector, blocks, adjoint=False):
    for v, block in blocks.items():
        vector = apply_block(vector, v, block, adjoint)
    return vector


def walsh_transform(vector):
    return [gsum(scale(value, character(x, y)) for y, value in enumerate(vector))
            for x in range(len(vector))]


def inner(a, b):
    return gsum(mul(x, conj(y)) for x, y in zip(a, b))


def physical_rows(order, blocks):
    values = []
    for x in range(order):
        value = ONE
        for v, block in blocks.items():
            value = mul(value, block if character(v, x) == 1 else conj(block))
        values.append(value)
    return values


def local_membership(vector, blocks):
    return all(divides(block, add(vector[y], vector[y ^ v]))
               and divides(conj(block), sub(vector[y], vector[y ^ v]))
               for v, block in blocks.items() for y in range(len(vector)))


def product_membership(vector, blocks):
    total_norm = prod(norm(block) for block in blocks.values())
    inverse_numerator = apply_all(vector, blocks, adjoint=True)
    return all(a % total_norm == b % total_norm == 0 for a, b in inverse_numerator)


def physical_membership(vector, blocks):
    transformed = walsh_transform(vector)
    return all(divides(row, value) for row, value in
               zip(physical_rows(len(vector), blocks), transformed))


def expansion(order, blocks):
    coefficients = [ZERO] * order
    labels = list(blocks)
    for choices in product((0, 1), repeat=len(labels)):
        label, value = 0, ONE
        for selected, v in zip(choices, labels):
            a, b = blocks[v]
            value = mul(value, (0, b) if selected else (a, 0))
            if selected:
                label ^= v
        coefficients[label] = add(coefficients[label], value)
    return coefficients


def core_fixtures():
    rng = Random(279)
    checked = 0
    for order in (4, 8, 16):
        primes = split_primes(order - 1)
        blocks = {v: power(primes[v - 1], 1 + (v % 2)) for v in range(1, order)}
        total_norm = prod(norm(block) for block in blocks.values())
        canonical = apply_all(delta(order), blocks)
        rows = physical_rows(order, blocks)
        assert walsh_transform(canonical) == rows
        assert walsh_transform(rows) == [scale(c, order) for c in canonical]
        if order <= 8:
            assert expansion(order, blocks) == canonical
        assert sum(map(norm, canonical)) == total_norm
        assert all(inner(translate(canonical, h), canonical) ==
                   ((total_norm, 0) if h == 0 else ZERO) for h in range(order))
        assert gsum(canonical) == rows[0]
        for bit in range(order.bit_length() - 1):
            assert scale(gsum(canonical[y] for y in range(order) if y & (1 << bit)), 2) == sub(rows[0], rows[1 << bit])
        common = (3, 2)  # May overlap a core prime; scalar content stays explicit.
        common_canonical = [mul(common, c) for c in canonical]
        source_norm = norm(common) * total_norm
        assert sum(map(norm, common_canonical)) == source_norm
        for h in range(order):
            basis_vector = [mul(common, c) for c in apply_all(delta(order, h), blocks)]
            assert basis_vector == translate(common_canonical, h)
            assert sum(map(norm, basis_vector)) == source_norm
        for _ in range(8):
            u = [(rng.randrange(-2, 3), rng.randrange(-2, 3)) for _ in range(order)]
            vector = apply_all(u, blocks)
            assert sum(map(norm, vector)) == total_norm * sum(map(norm, u))
            assert apply_all(vector, blocks, adjoint=True) == [scale(value, total_norm) for value in u]
            for perturbation in (ZERO, ONE, (0, 1)):
                candidate = vector[:]
                candidate[0] = add(candidate[0], perturbation)
                assert local_membership(candidate, blocks) == product_membership(candidate, blocks)
                assert physical_membership(candidate, blocks) == product_membership(candidate, blocks)
                checked += 1
    return checked


def single_block_exhaustion():
    checked = 0
    block = (2, 1)
    for a, b, c, d in product(range(-3, 4), repeat=4):
        vector = [(a, b), (c, d)]
        assert local_membership(vector, {1: block}) == product_membership(vector, {1: block})
        checked += 1
    return checked


def overlap_counterexample():
    blocks = {1: (2, 1), 2: (3, 2), 3: (2, 1)}
    vector = apply_block(delta(4, value=(5, 0)), 2, blocks[2])
    total_norm = prod(norm(block) for block in blocks.values())
    assert sum(map(norm, vector)) == total_norm
    assert local_membership(vector, blocks)
    assert not product_membership(vector, blocks)
    assert not physical_membership(vector, blocks)
    inverse = [scale(value, Fraction(1, total_norm))
               for value in apply_all(vector, blocks, adjoint=True)]
    assert inverse == [(Fraction(4, 5), 0), (0, Fraction(-2, 5)),
                       (Fraction(-1, 5), 0), (0, Fraction(-2, 5))]
    source_labels = (0, 1, 2)
    shared_allocations = [sum(character(v, x) == 1 for v in (1, 3))
                          for x in source_labels]
    assert shared_allocations == [2, 0, 1]
    assert {x for x, a in zip(source_labels, shared_allocations) if a >= 1} == {0, 2}
    assert {x for x, a in zip(source_labels, shared_allocations) if a >= 2} == {0}
    # Omitting a ramified common-factor condition also gives a false test.
    assert all(divides((2, 0), z) for z in walsh_transform([ONE] * 4))
    assert not all(divides((2, 0), z) for z in [ONE] * 4)


def unit_patterns_and_corrections():
    blocks = {1: (2, 1), 2: (3, 2), 3: (4, 1)}
    canonical = apply_all(delta(4), blocks)
    rows = physical_rows(4, blocks)
    accepted = 0
    for units in product(UNITS, repeat=4):
        is_character = any(all(units[x] == scale(unit, character(h, x)) for x in range(4))
                           for unit in UNITS for h in range(4))
        transformed_units = walsh_transform(units)
        integral_operator = all(a % 4 == b % 4 == 0 for a, b in transformed_units)
        assert integral_operator == is_character
        transformed_rows = walsh_transform([mul(unit, row) for unit, row in zip(units, rows)])
        integral_coefficients = all(a % 4 == b % 4 == 0 for a, b in transformed_rows)
        assert integral_coefficients == is_character  # Odd core determinant.
        if integral_coefficients:
            vector = [(a // 4, b // 4) for a, b in transformed_rows]
            assert product_membership(vector, blocks)
            assert any(vector == [mul(unit, c) for c in translate(canonical, h)]
                       for unit in UNITS for h in range(4))
            accepted += 1
    assert accepted == 16
    corrections = [(2, 1)] * 4
    corrections[1] = (2, -1)
    assert all(norm(c) == 5 for c in corrections)
    transformed = walsh_transform(corrections)
    assert not all(a % 4 == b % 4 == 0 for a, b in transformed)
    units = [ONE, (0, 1), (-1, 0), (0, -1)]
    actual = [mul(mul(unit, correction), row)
              for unit, correction, row in zip(units, corrections, rows)]
    for bit in (0, 1):
        index = 1 << bit
        cut_sum = gsum(canonical[y] for y in range(4) if y & index)
        first_correction = mul(units[0], corrections[0])
        second_correction = mul(units[index], corrections[index])
        center = scale(mul(actual[0], sub(gdiv(ONE, first_correction),
                                          gdiv(ONE, second_correction))), Fraction(1, 2))
        error = sub(cut_sum, center)
        assert error == gdiv(sub(actual[0], actual[index]), scale(second_correction, 2))
        assert norm(error) == Fraction(norm(sub(actual[0], actual[index])), 20)
    return accepted


if __name__ == "__main__":
    memberships = core_fixtures()
    local = single_block_exhaustion()
    overlap_counterexample()
    units = unit_patterns_and_corrections()
    print("PASS: literal order-4, 8, and 16 Fourier coefficients, energy, autocorrelation, and orthogonal bases.")
    print(f"PASS: {memberships} full membership checks and {local} exhaustive single-block congruence checks.")
    print("PASS: overlapping-block false minimum and the separate ramified common-content condition.")
    print(f"PASS: all 256 order-4 unit functions; exactly {units} unit/character twists preserve integrality.")
    print("PASS: equal-norm correction counterexample and exact source-error phase centers.")
    print("No endpoint example, moving-core height gain, or universal phase-propagation claim is made.")
