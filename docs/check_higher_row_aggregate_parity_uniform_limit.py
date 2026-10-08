"""Exact coefficient and parity audit for the higher-row aggregate argument."""

from fractions import Fraction
from itertools import combinations, product
from math import comb


def balanced_representatives(k):
    vectors = [v for v in product((-1, 1), repeat=k) if sum(v) == 0]
    return [v for v in vectors if v[0] == 1]


def coefficient_from_characters(k, r):
    chars = balanced_representatives(k)
    # The canonical cut has r minus signs.  Absolute values make its
    # orientation and the choice of representative immaterial.
    cut = (-1,) * r + (1,) * (k - r)
    return sum(Fraction(abs(sum(x * y for x, y in zip(lam, cut))), 4)
               for lam in chars) / len(chars)


def coefficient_hypergeometric(k, r):
    total = comb(k, k // 2)
    numerator = 0
    for a in range(r + 1):
        if 0 <= k // 2 - a <= k - r:
            numerator += (comb(r, a) * comb(k - r, k // 2 - a)
                          * abs(r - 2 * a))
    return Fraction(numerator, 2 * total)


for k in (4, 6, 8, 10, 12):
    chars = balanced_representatives(k)
    assert len(chars) == comb(k, k // 2) // 2
    # Since every lambda_i is odd modulo two, signed-exponent parity is
    # independent of the balanced character.
    for parity_vector in product((0, 1), repeat=k):
        expected_parity = sum(parity_vector) % 2
        assert all(sum(x * y for x, y in zip(lam, parity_vector)) % 2
                   == expected_parity for lam in chars)
    for r in range(1, k // 2 + 1):
        assert coefficient_from_characters(k, r) == coefficient_hypergeometric(k, r)

expected = {
    4: (Fraction(1, 2), Fraction(1, 3)),
    6: (Fraction(1, 2), Fraction(2, 5), Fraction(3, 5)),
    8: (Fraction(1, 2), Fraction(3, 7), Fraction(9, 14), Fraction(18, 35)),
    10: (Fraction(1, 2), Fraction(4, 9), Fraction(2, 3), Fraction(4, 7), Fraction(5, 7)),
}
for k, values in expected.items():
    assert tuple(coefficient_hypergeometric(k, r)
                 for r in range(1, k // 2 + 1)) == values

# The six-row profile has 6, 15, 10 unoriented cuts at sizes 1, 2, 3.
k = 6
counts = [comb(k, r) if r < k // 2 else comb(k, r) // 2
          for r in range(1, k // 2 + 1)]
coefficients = [coefficient_hypergeometric(k, r)
                for r in range(1, k // 2 + 1)]
assert counts == [6, 15, 10]
assert sum(n * c for n, c in zip(counts, coefficients)) == 15
assert sum(counts) == 31
assert sum(n * (c - Fraction(1, 2))
           for n, c in zip(counts, coefficients)) == Fraction(-1, 2)

# A 1/20 relative cut-weight neighborhood remains strictly negative.
eta = Fraction(1, 20)
upper = (Fraction(1, 10) *
         (10 * (1 + eta) - 15 * (1 - eta)))
assert upper == Fraction(-3, 8)
assert Fraction(6, 8) == Fraction(3, 4)  # six-row angular factor k/8

# Explicit unbounded-support parity counterfamily: row i is odd only at
# its private prime p_i.  No six-subset has zero coordinatewise parity.
for m in range(6, 25):
    rows = [1 << i for i in range(m)]
    assert all(sum(rows[i] for i in subset) != 0
               for subset in combinations(range(m), 6))


def canonical_cut(mask, m):
    full = (1 << m) - 1
    complement = full ^ mask
    return min(mask, complement)


def integer_rank(rows):
    rows = [[Fraction(x) for x in row] for row in rows]
    rank = 0
    columns = len(rows[0]) if rows else 0
    for column in range(columns):
        pivot = next((i for i in range(rank, len(rows))
                      if rows[i][column]), None)
        if pivot is None:
            continue
        rows[rank], rows[pivot] = rows[pivot], rows[rank]
        pivot_value = rows[rank][column]
        rows[rank] = [x / pivot_value for x in rows[rank]]
        for i in range(len(rows)):
            if i == rank or rows[i][column] == 0:
                continue
            factor = rows[i][column]
            rows[i] = [x - factor*y for x, y in zip(rows[i], rows[rank])]
        rank += 1
        if rank == len(rows):
            break
    return rank


# Add dominant even-exponent columns for every nonconstant cut, together
# with one private odd column per row.  This keeps all cut weights nearly
# uniform while preventing every six-row parity sum from vanishing.
for m in (6, 8):
    full = (1 << m) - 1
    cuts = [mask for mask in range(1, full) if mask < (full ^ mask)]
    for q in (1, 10, 100):
        allocations = []
        cut_weights = {mask: 2*q for mask in cuts}
        for mask in cuts:
            allocations.append(tuple(2*q if mask >> i & 1 else 0
                                     for i in range(m)))
        for i in range(m):
            allocations.append(tuple(1 if j == i else 0 for j in range(m)))
            cut_weights[canonical_cut(1 << i, m)] += 1

        assert max(cut_weights.values()) - min(cut_weights.values()) <= 1
        assert all(abs(cut_weights[mask] - 2*q) <= 1 for mask in cuts)
        assert integer_rank(allocations) == m
        total_weight = sum(cut_weights.values())
        for i, j in combinations(range(m), 2):
            assert any(a[i] != a[j] for a in allocations)
            distance = sum(weight for cut, weight in cut_weights.items()
                           if ((cut >> i) ^ (cut >> j)) & 1)
            assert Fraction(distance) - Fraction(total_weight, 2) == q + 2 - Fraction(m, 2)
            if q > Fraction(m - 4, 2):
                assert 2 * distance > total_weight
        if m >= 6:
            parity = [1 << i for i in range(m)]
            assert all(sum(parity[i] for i in subset) != 0
                       for subset in combinations(range(m), 6))

print("PASS: exact k=4,6,8,10,12 coefficient tables; six-row margin -1/62;")
print("      1/20-neighborhood bound -3/8; unbounded-support parity obstruction;")
print("      dominant even full-cut columns preserve near-uniform profiles and rank.")
