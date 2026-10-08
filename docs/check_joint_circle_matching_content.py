"""Exact local cancellation example for the joint Gaussian circle lift.

This is an actual integral six-tuple at one split prime. It is not a
full-profile endpoint example or a lattice-point counterexample.
"""

from itertools import combinations, product
from math import gcd


def add(a, b):
    return a[0] + b[0], a[1] + b[1]


def scale(a, n):
    return a[0] * n, a[1] * n


def mul(a, b):
    return a[0] * b[0] - a[1] * b[1], a[0] * b[1] + a[1] * b[0]


def conjugate(a):
    return a[0], -a[1]


def norm(a):
    return a[0] ** 2 + a[1] ** 2


def exact_divide(a, b):
    numerator = mul(a, conjugate(b))
    denominator = norm(b)
    assert all(x % denominator == 0 for x in numerator)
    return tuple(x // denominator for x in numerator)


def gaussian_gcd(a, b):
    while b != (0, 0):
        numerator = mul(a, conjugate(b))
        denominator = norm(b)
        quotient = tuple((2 * x + denominator) // (2 * denominator)
                         for x in numerator)
        remainder = add(a, scale(mul(quotient, b), -1))
        assert norm(remainder) < norm(b)
        a, b = b, remainder
    return a


def valuation(a, pi):
    assert a != (0, 0)
    answer = 0
    while True:
        numerator = mul(a, conjugate(pi))
        if any(x % norm(pi) for x in numerator):
            return answer
        a = tuple(x // norm(pi) for x in numerator)
        answer += 1


def determinant(a, b):
    return a[0] * b[1] - a[1] * b[0]


def matchings(indices):
    if not indices:
        yield ()
        return
    first = indices[0]
    for second in indices[1:]:
        for tail in matchings(tuple(i for i in indices if i not in (first, second))):
            yield ((first, second),) + tail


MATCHINGS = tuple(matchings(tuple(range(6))))
assert len(MATCHINGS) == 15

# The sorted-valuation formula has a proof in the accompanying note.
for orders in product(range(4), repeat=6):
    lower = min(sum(min(orders[i], orders[j]) for i, j in matching)
                for matching in MATCHINGS)
    assert lower == sum(sorted(orders)[:3])

P = [(x, 1) for x in (8, 2, 0, 21, 1, 9)]


def d(i, j):
    return determinant(P[i - 1], P[j - 1])


A = d(1, 2) * d(3, 4) * d(5, 6)
B = d(1, 2) * d(3, 6) * d(4, 5)
C = d(1, 4) * d(2, 3) * d(5, 6)
F, H, V = C - B, C + B, 2 * A + B + C
assert (A, B, C, F, H, V) == (1008, -1080, 208, 1288, -872, 1144)

Q = [add(scale(P[3], H), scale(P[1], -2 * d(3, 4) * d(5, 6) * d(4, 1))),
     add(scale(P[4], F), scale(P[1], -2 * d(3, 5) * d(4, 6) * d(5, 1))),
     add(scale(P[5], H), scale(P[1], 2 * d(3, 6) * d(4, 5) * d(6, 1)))]
assert Q == [(-27048, -5240), (952, 1120), (-8568, -1232)]

pi = (3, 2)
assert [valuation(p, pi) for p in P] == [1, 0, 0, 1, 0, 0]
assert [valuation(p, conjugate(pi)) for p in P] == [0] * 6
assert [(i, j) for i, j in combinations(range(1, 7), 2) if d(i, j) % 13 == 0] == [(1, 4)]
assert valuation((d(1, 4), 0), pi) == 1
assert [valuation(q, pi) for q in Q] == [1, 1, 0]
assert [valuation(q, conjugate(pi)) for q in Q] == [0, 0, 0]

W = P[:3] + Q
determinant_products = []
for matching in MATCHINGS:
    value = 1
    for i, j in matching:
        value *= determinant(W[i], W[j])
    determinant_products.append(value)
ordinary_content = 0
for value in determinant_products:
    ordinary_content = gcd(ordinary_content, value)
assert ordinary_content % 13 != 0

# Build an actual equal-norm integral circle, retaining anchor z0.
fractions = []
common_denominator = (1, 0)
for w in W:
    content = gaussian_gcd(w, conjugate(w))
    numerator = exact_divide(w, content)
    denominator = exact_divide(conjugate(w), content)
    fractions.append((numerator, denominator))
    overlap = gaussian_gcd(common_denominator, denominator)
    common_denominator = mul(exact_divide(common_denominator, overlap), denominator)
points = [mul(exact_divide(common_denominator, denominator), numerator)
          for numerator, denominator in fractions]
assert len(set(points + [common_denominator])) == 7
assert all(norm(point) == norm(common_denominator) for point in points)

gamma = (0, 0)
for matching in MATCHINGS:
    value = (1, 0)
    for i, j in matching:
        value = mul(value, add(points[i], scale(points[j], -1)))
    gamma = gaussian_gcd(gamma, value)
assert valuation(gamma, pi) == valuation(gamma, conjugate(pi)) == 0

print("Sorted-valuation matching identity checked on all4096 test tuples.")
print("Exact S14 example: old phase cut(1,4) becomes balanced output cut(1,4,5).")
print("Actual seven-point equal-norm integral realization checked.")
print("Its common matching scalar has valuation zero at both primes over13.")
print("No global full-profile realization or uniform-bound counterexample is asserted.")
