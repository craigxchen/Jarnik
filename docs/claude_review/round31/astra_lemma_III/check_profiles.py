"""Exact finite checks for the five-mark proofs, independent of Gaussian angles."""

from fractions import Fraction
from itertools import combinations, product
from math import gcd

from sympy import Matrix


marks = tuple(range(5))
pairs = tuple(combinations(marks, 2))
cuts = [frozenset([i]) for i in marks] + list(map(frozenset, pairs))


def crosses(cut, pair):
    return int((pair[0] in cut) != (pair[1] in cut))


def anchored_column(cut):
    return Matrix([int(i in cut) - int(4 in cut) for i in range(4)])


def column_index(columns):
    determinants = (abs(int(Matrix.hstack(*subset).det()))
                    for subset in combinations(columns, 4))
    return gcd(*determinants)


assert column_index([anchored_column(e) for e in pairs]) == 1

for i, j in pairs:
    tested = [e for e in pairs if i in e or j in e]
    assert len(tested) == 7
    for cut in cuts:
        count = sum(crosses(cut, e) for e in tested)
        if cut == frozenset([i, j]):
            assert count == 6
        else:
            assert count <= 4

# Any missing disjoint pair of edges gives a four-cycle certificate.
for e, f in combinations(pairs, 2):
    if set(e).isdisjoint(f):
        cycle = tuple(product(e, f))
        assert all(sum(crosses(cut, x) for x in cycle) <= 2 for cut in cuts
                   if cut not in (set(e), set(f)))

def word_cost(weights):
    return sum(abs(weights[i] + weights[j]) for i, j in pairs)


assert word_cost((1, -1, 0, 0, 0)) == 6
for pattern, cost in [((2, -2, 0, 0, 0), 12),
                      ((2, -1, -1, 0, 0), 12),
                      ((1, 1, -1, -1, 0), 8)]:
    assert word_cost(pattern) == cost

quartet_matchings = [((0, 1), (2, 3)), ((0, 2), (1, 3)), ((0, 3), (1, 2))]
for a in product((0, 1), repeat=4):
    cost = sum(abs(a[i] + a[j] - a[k] - a[l])
               for (i, j), (k, l) in quartet_matchings)
    assert cost == (0, 3, 2, 3, 0)[sum(a)]

for word, expected in [((1, -1, 0, 0, 0), 8), ((1, 1, -1, -1, 0), 12)]:
    assert sum(abs(sum(word[i] for i in cut)) for cut in cuts) == expected

# Comparisons for the volume bound, with all radicals removed.
assert 70 * 5 > (Fraction(5, 4) * Fraction(17, 5)) ** 4
assert 70 * 7 > (Fraction(5, 4) * Fraction(94, 25)) ** 4
assert Fraction(94, 25) > Fraction(124, 33) > Fraction(15, 4)
print("Pair-cut index 1, cut certificates, square-word costs and thresholds verified.")

# Barycentric identity (8) on a genuine ten-prime tuple.
blocks = [(2, 1), (3, 2), (4, 1), (5, 2), (6, 1),
          (5, 4), (7, 2), (6, 5), (8, 3), (8, 5)]
norms = [a*a + b*b for a, b in blocks]
assert norms == [5, 13, 17, 29, 37, 41, 53, 61, 73, 89]
points = [(1, 0) for _ in marks]
for (i, j), (a, b) in zip(pairs, blocks):
    for k in (i, j):
        x, y = points[k]
        points[k] = (x*a - y*b, x*b + y*a)
denominators = [1 for _ in marks]
for (i, j), n in zip(pairs, norms):
    x, y = points[i]
    u, v = points[j]
    quotient, remainder = divmod(x*v - y*u, n)
    assert remainder == 0 and quotient != 0
    denominators[i] *= quotient
    denominators[j] *= -quotient
for coordinate in (0, 1):
    assert sum(Fraction(point[coordinate], d)
               for point, d in zip(points, denominators)) == 0
print("Both coordinates of the additive identity verified on a genuine ten-prime tuple.")
