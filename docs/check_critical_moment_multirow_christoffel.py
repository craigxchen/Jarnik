#!/usr/bin/env python3
"""Exact algebra audit, not a distinct-node odd-moment solution.

The partial Walsh fixture verifies the centered-row projector. Arbitrary
distinct positive square nodes verify the Hankel projector, the left-side
identity in (6), and Cauchy--Binet. These nodes do NOT solve CU=0; consequently
this checker does not assert K+Q<=I or inequality (6) on that fixture.
Only Python's standard library is required.
"""

from fractions import Fraction as F
from itertools import combinations, product


def transpose(a):
    return [list(row) for row in zip(*a)]


def multiply(a, b):
    return [[sum(x * y for x, y in zip(row, col))
             for col in transpose(b)] for row in a]


def determinant(a):
    a = [[F(x) for x in row] for row in a]
    result = F(1)
    for j in range(len(a)):
        pivot = next((i for i in range(j, len(a)) if a[i][j]), None)
        if pivot is None:
            return F(0)
        if pivot != j:
            a[j], a[pivot] = a[pivot], a[j]
            result = -result
        value = a[j][j]
        result *= value
        for i in range(j + 1, len(a)):
            scale = a[i][j] / value
            a[i] = [x - scale * y for x, y in zip(a[i], a[j])]
    return result


def inverse(a):
    n = len(a)
    a = [[F(x) for x in row] + [F(i == j) for j in range(n)]
         for i, row in enumerate(a)]
    for j in range(n):
        pivot = next(i for i in range(j, n) if a[i][j])
        a[j], a[pivot] = a[pivot], a[j]
        value = a[j][j]
        a[j] = [x / value for x in a[j]]
        for i in range(n):
            if i != j:
                scale = a[i][j]
                a[i] = [x - scale * y for x, y in zip(a[i], a[j])]
    return [row[n:] for row in a]


def prod(values):
    result = 1
    for value in values:
        result *= value
    return result


def main():
    s, n = 3, 14
    t = [[(-1) ** bin(i & j).count('1') for j in range(16)]
         for i in range(8)]
    assert multiply(t, transpose(t)) == [
        [16 * (i == j) for j in range(8)] for i in range(8)]
    b = [row[1] for row in t]
    v = [row[2:] for row in t]
    h = [sum(row[j] for row in v) for j in range(n)]
    d = [sum(b[i] * v[i][j] for i in range(8)) for j in range(n)]
    c = [[v[i][j] - F(h[j], 8) for j in range(n)] for i in range(8)]
    q = [[F(sum(v[i][j] * v[i][k] for i in range(8)), 4 * (s + 1))
          - F(h[j] * h[k], 32 * (s + 1))
          + F(d[j] * d[k], 16 * (s * s - 1))
          for k in range(n)] for j in range(n)]
    assert q == transpose(q)
    assert multiply(q, q) == q
    assert sum(q[j][j] for j in range(n)) == 7
    assert multiply(c, q) == c
    assert multiply(c, transpose(c)) == [
        [16 * (F(i == k) - F(1, 8)) - b[i] * b[k]
         for k in range(8)] for i in range(8)]
    for j in range(n):
        ha = sum(v[i][j] for i in range(8) if b[i] == 1)
        hb = sum(v[i][j] for i in range(8) if b[i] == -1)
        diagonal = (8 - F(ha * ha + hb * hb, 4)) / (4 * (s + 1))
        diagonal += F((ha - hb) ** 2, 32 * (s - 1))
        assert q[j][j] == diagonal

    # Positive square nodes keep U and both projectors rational.
    roots = list(range(1, n + 1))
    x = [r * r for r in roots]
    u = [[roots[j] * x[j] ** k for k in range(s)] for j in range(n)]
    hankel = [[sum(z ** (k + ell + 1) for z in x)
               for ell in range(s)] for k in range(s)]
    assert multiply(transpose(u), u) == hankel
    assert all(determinant([row[:m] for row in hankel[:m]]) > 0
               for m in range(1, s + 1))
    kernel = multiply(multiply(u, inverse(hankel)), transpose(u))
    assert kernel == transpose(kernel)
    assert multiply(kernel, kernel) == kernel
    assert sum(kernel[j][j] for j in range(n)) == s
    assert multiply(kernel, u) == u
    assert any(value for row in multiply(c, u) for value in row)
    det_h = determinant(hankel)
    total = 0
    subset_count = 0
    for subset in combinations(range(n), s):
        weight = prod(x[j] for j in subset)
        weight *= prod((x[k] - x[j]) ** 2 for j, k in combinations(subset, 2))
        minor = [[kernel[j][k] for k in subset] for j in subset]
        assert determinant(minor) == F(weight, det_h)
        total += weight
        subset_count += 1
    assert total == det_h

    # For all s>=2, multiply q-q_min by 32(s^2-1). Its numerator
    # is affine in s. Check its slope and its value at s=2 exactly.
    checked_types = 0
    for ha, hb in product(range(-4, 5, 2), repeat=2):
        if (ha, hb) in [(-4, -4), (4, 4)]:
            continue
        slope = 36 - (ha + hb) ** 2
        constant = -44 + 3 * ha * ha + 3 * hb * hb - 2 * ha * hb
        assert slope >= 0 and 2 * slope + constant >= 0
        checked_types += 1
    assert checked_types == 23
    # q_min - 1/(2s+1) has positive numerator 6s^2-3s+3.
    # Writing s=2+t gives 6t^2+21t+21, all coefficients positive.
    assert (6, 24 - 3, 24 - 6 + 3) == (6, 21, 21)
    print('PASS: rank-seven Walsh projector; rational Hankel projector;')
    print('      {} Vandermonde minors and Cauchy--Binet;'.format(subset_count))
    print('      all 23 nonconstant group-sum types, uniformly for s>=2.')
    print('SCOPE: arbitrary nodes do not solve CU=0; no moment-solution or PSD feasibility claim.')


if __name__ == '__main__':
    main()
