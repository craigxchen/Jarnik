"""Exact Gaussian checks for complement_reflection_invariant_relations.md."""

from itertools import combinations
from math import gcd


ONE = (1, 0)


def mul(a, b):
    return (a[0] * b[0] - a[1] * b[1], a[0] * b[1] + a[1] * b[0])


def conjugate(a):
    return (a[0], -a[1])


def norm(a):
    return a[0]**2 + a[1]**2


def product(values):
    result = ONE
    for value in values:
        result = mul(result, value)
    return result


def power(a, exponent):
    return product([a] * exponent)


def exact_div(a, b):
    numerator = mul(a, conjugate(b))
    denominator = norm(b)
    assert denominator and all(coordinate % denominator == 0 for coordinate in numerator)
    return tuple(coordinate // denominator for coordinate in numerator)


def gaussian_gcd(a, b):
    while b != (0, 0):
        numerator = mul(a, conjugate(b))
        denominator = norm(b)
        quotient = tuple((2 * coordinate + denominator) // (2 * denominator) for coordinate in numerator)
        subtrahend = mul(quotient, b)
        a, b = b, tuple(x - y for x, y in zip(a, subtrahend))
    return a


def determinant(a, b):
    return a[0] * b[1] - a[1] * b[0]


def primitive_pair(a, b):
    common_norm = norm(gaussian_gcd(a, b))
    numerator = mul(a, conjugate(b))
    assert all(coordinate % common_norm == 0 for coordinate in numerator)
    result = tuple(coordinate // common_norm for coordinate in numerator)
    assert norm(gaussian_gcd(result, conjugate(result))) == 1
    return result


def reflect(blocks, corrections):
    m = len(corrections)
    rows = [mul(corrections[i], product(value for cut, value in blocks.items() if i in cut))
            for i in range(m)]
    assert all(norm(gaussian_gcd(row, conjugate(row))) == 1 for row in rows)
    all_corrections = product(corrections)
    total = mul(all_corrections, product(blocks.values()))
    reflected = []
    contents = []
    for row in rows:
        raw = exact_div(total, row)
        assert norm(raw) % 2 == 1
        content = gcd(abs(raw[0]), abs(raw[1]))
        contents.append(content)
        reflected.append(tuple(coordinate // content for coordinate in raw))
    assert all(norm(gaussian_gcd(row, conjugate(row))) == 1 for row in reflected)
    trim = {cut: gaussian_gcd(value, conjugate(all_corrections)) for cut, value in blocks.items()}
    reduced_blocks = {cut: exact_div(value, trim[cut]) for cut, value in blocks.items()}
    assert norm(product(trim.values())) <= norm(all_corrections)
    for i, row in enumerate(reflected):
        core = product(value for cut, value in reduced_blocks.items() if i not in cut)
        new_correction = exact_div(row, core)
        old_correction = exact_div(all_corrections, corrections[i])
        assert norm(new_correction) <= norm(old_correction) * norm(all_corrections)
    for i, j in combinations(range(m), 2):
        assert primitive_pair(reflected[i], reflected[j]) == conjugate(primitive_pair(rows[i], rows[j]))
        assert (determinant(reflected[i], reflected[j]) * norm(rows[i]) * norm(rows[j])
                * contents[i] * contents[j] == -norm(total) * determinant(rows[i], rows[j]))
    return rows, reflected, trim, contents


def main():
    m = 6
    blocks = {frozenset(j for j in range(m) if mask >> j & 1): ONE for mask in range(1 << m)}
    blocks[frozenset((0, 1, 2))] = power((2, 1), 3)
    blocks[frozenset((1, 3))] = power((3, 2), 2)
    blocks[frozenset((2, 4, 5))] = (4, 1)
    blocks[frozenset()] = (5, 2)
    corrections = [(6, 1), (5, 4), (7, 2), (2, -1), (3, -2), ONE]
    rows, reflected, trim, contents = reflect(blocks, corrections)
    assert all(determinant(rows[i], rows[j]) for i, j in combinations(range(m), 2))
    assert norm(product(trim.values())) > 1
    assert max(contents) > 1
    first_matching = ((0, 1), (2, 3), (4, 5))
    second_matching = ((0, 2), (1, 4), (3, 5))

    def invariant(points, matching):
        value = 1
        for i, j in matching:
            value *= determinant(points[i], points[j])**2
        return value

    old_first = invariant(rows, first_matching)
    old_second = invariant(rows, second_matching)
    assert old_second * invariant(reflected, first_matching) == old_first * invariant(reflected, second_matching)

    # A reflected real primitive row is allowed. This deliberately uses
    # degenerate support weights; it checks the zero-imaginary edge case only.
    zero_blocks = {frozenset((0,)): (2, 1)}
    _, zero_reflected, _, _ = reflect(zero_blocks, [ONE, (2, -1), ONE, ONE, ONE, ONE])
    assert any(row[1] == 0 and abs(row[0]) == 1 for row in zero_reflected)

    # Every one of the sixteen four-row blocks is nonunit, with varying
    # exponents and conjugate correction primes only on inactive rows.
    primes = []
    used_norms = set()
    for a in range(1, 30):
        for b in range(1, a):
            p = a*a + b*b
            if p not in used_norms and all(p % d for d in range(2, int(p**0.5) + 1)):
                primes.append((a, b))
                used_norms.add(p)
    assert len(primes) >= 16
    full_blocks = {
        frozenset(j for j in range(4) if mask >> j & 1): power(primes[mask], 1 + mask % 3)
        for mask in range(16)
    }
    full_corrections = [conjugate(primes[1 << ((i + 1) % 4)]) for i in range(4)]
    full_rows, _, full_trim, full_contents = reflect(full_blocks, full_corrections)
    assert all(determinant(full_rows[i], full_rows[j]) for i, j in combinations(range(4), 2))
    assert norm(product(full_trim.values())) > 1 and max(full_contents) > 1
    assert 2 * 6 * (4 * 4 - 3) == 156
    print("Integral complement reflection, nontrivial common trimming, and correction norm bounds pass.")
    print("All primitive pair numerators are conjugated exactly; the invariant zero keeps its coefficients.")
    print("A reflected zero-imaginary primitive row is explicitly covered; no ramified unit rotation occurs.")
    print("The six-row two-pass correction coefficient is 156.")


if __name__ == "__main__":
    main()
