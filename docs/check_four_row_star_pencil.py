#!/usr/bin/env python3
"""Exact interval audit of the star-pencil form of the certified root.

The existence/nondegeneracy proof is the separate contraction certificate.
This checker uses its rational coordinate cube, never a reconstructed root.
"""
from fractions import Fraction as F
from itertools import permutations
import json
from pathlib import Path


class I:
    def __init__(self, lo=0, hi=None):
        self.lo = F(lo)
        self.hi = self.lo if hi is None else F(hi)
        assert self.lo <= self.hi

    def __add__(self, other):
        other = other if isinstance(other, I) else I(other)
        return I(self.lo + other.lo, self.hi + other.hi)

    __radd__ = __add__

    def __neg__(self):
        return I(-self.hi, -self.lo)

    def __sub__(self, other):
        return self + (-other if isinstance(other, I) else -I(other))

    def __mul__(self, other):
        other = other if isinstance(other, I) else I(other)
        candidates = [a * b for a in (self.lo, self.hi)
                      for b in (other.lo, other.hi)]
        return I(min(candidates), max(candidates))

    __rmul__ = __mul__

    def __truediv__(self, other):
        other = other if isinstance(other, I) else I(other)
        assert other.lo > 0 or other.hi < 0
        return self * I(1 / other.hi, 1 / other.lo)

    def nonzero(self):
        return self.lo > 0 or self.hi < 0

    def midpoint(self):
        return (self.lo + self.hi) / 2


def cmul(a, b):
    return a[0] * b[0] - a[1] * b[1], a[0] * b[1] + a[1] * b[0]


def conj(a):
    return a[0], -a[1]


def polynomial(roots):
    coefficients = [(I(1), I(0))]
    for root in roots:
        nxt = [(I(0), I(0)) for _ in range(len(coefficients) + 1)]
        for k, coefficient in enumerate(coefficients):
            product = cmul(coefficient, root)
            nxt[k] = tuple(a - b for a, b in zip(nxt[k], product))
            nxt[k + 1] = tuple(a + b for a, b in zip(nxt[k + 1], coefficient))
        coefficients = nxt
    return coefficients


def determinant(matrix):
    result = I(0)
    for permutation in permutations(range(len(matrix))):
        inversions = sum(permutation[a] > permutation[b]
                         for a in range(len(matrix))
                         for b in range(a + 1, len(matrix)))
        term = I((-1) ** inversions)
        for a, b in enumerate(permutation):
            term = term * matrix[a][b]
        result = result + term
    return result


def main():
    data = json.loads(Path(__file__).with_name(
        "four_row_real_polynomial_certificate.json").read_text())
    rho = F(data["radius"])
    assert rho == F(1, 10**30)
    roots = [tuple(I(F(c) - rho, F(c) + rho) for c in root)
             for root in data["roots"]]
    roots[14] = I(0), I(1)
    vertices = [conj(roots[14])] + [roots[(1 << i) - 1] for i in range(4)]
    edges = {(a, b): conj(roots[(15 ^ (1 << (b - 1))) - 1]) if a == 0
             else roots[((1 << (a - 1)) | (1 << (b - 1))) - 1]
             for a in range(5) for b in range(a + 1, 5)}
    stars = [polynomial([vertices[a]] +
                        [edges[min(a, b), max(a, b)] for b in range(5) if b != a])
             for a in range(5)]
    coefficient_rows = [[stars[a][k][0] for a in range(5)]
                        for k in (5, 4, 3, 2)]
    cofactors = [(-1) ** a * determinant(
        [[row[b] for b in range(5) if b != a] for row in coefficient_rows])
        for a in range(5)]
    assert all(value.nonzero() for value in cofactors)
    assert all(value.lo > 0 if sign > 0 else value.hi < 0
               for value, sign in zip(cofactors, (-1, 1, 1, -1, -1)))
    assert cofactors[4].hi < -F(6, 1000)

    # M = [Re Q; Im Q].  Its degree-2 minors imply proportional leading
    # coefficient vectors whenever the two row degrees sum to more than 2.
    # Therefore each forced zero below is an exact symbolic consequence,
    # not a numerical truncation.  The intervals certify every pivot.
    matrix = [[[stars[a][k][part] for a in range(5)] for k in range(6)]
              for part in range(2)]
    matrix[1][5] = [I(0) for _ in range(5)]  # monic stars
    transform = [[[I(int(row == col and k == 0)) for col in range(2)]
                  for k in range(6)] for row in range(2)]
    degrees = [5, 4]
    degree_history = [tuple(degrees)]
    while sum(degrees) > 2:
        assert all(matrix[row][degrees[row]][0].nonzero() for row in range(2))
        row = int(degrees[0] < degrees[1])
        other = 1 - row
        shift = degrees[row] - degrees[other]
        factor = matrix[row][degrees[row]][0] / matrix[other][degrees[other]][0]
        for k in range(shift, 6):
            for a in range(5):
                matrix[row][k][a] = matrix[row][k][a] - factor * matrix[other][k - shift][a]
            for a in range(2):
                transform[row][k][a] = transform[row][k][a] - factor * transform[other][k - shift][a]
        matrix[row][degrees[row]] = [I(0) for _ in range(5)]
        degrees[row] -= 1
        degree_history.append(tuple(degrees))
    assert degrees == [1, 1]
    assert all(matrix[row][1][0].nonzero() for row in range(2))
    # U^{-1} gives Q=h*A+g*B.  Degree upper bounds follow directly from
    # these seven elementary row operations, and these intervals give equality.
    assert transform[1][4][1].nonzero()  # Re h, degree 4
    assert transform[0][3][1].nonzero()  # Re g, degree 3
    for row, limit in ((0, 3), (1, 4)):
        assert all(value.lo == value.hi == 0
                   for k in range(limit + 1, 6) for value in transform[row][k])
    assert transform[1][4][0].lo == transform[1][4][0].hi == 0
    assert transform[0][3][0].lo == transform[0][3][0].hi == 0
    common_slope = matrix[0][1][0]
    alphas = [matrix[0][0][a] / common_slope for a in range(5)]
    ds = [common_slope * matrix[1][0][a]
          - matrix[0][0][a] * matrix[1][1][a] for a in range(5)]
    assert all(d.lo > 0 for d in ds)
    increasing_y = [0, 3, 1, 2, 4]
    assert all(alphas[a].lo > alphas[b].hi
               for a, b in zip(increasing_y, increasing_y[1:]))
    y_constants = [common_slope * (matrix[1][1][a] - matrix[1][1][0])
                   for a in range(5)]
    norm_pair_matrix = [[((y_constants[b] - y_constants[a])
                         * (alphas[b] - alphas[a]) + ds[a] + ds[b]) / 2
                        for b in range(5)] for a in range(5)]
    principal_minor = determinant([row[:3] for row in norm_pair_matrix[:3]])
    assert principal_minor.hi < -F(1, 10**9)
    print("PASS: certified star span has real rank exactly 4")
    print("relation cofactors:", [float(x.midpoint()) for x in cofactors])
    print("certified row degree sequence:", degree_history)
    print("PASS: Q_a=h*A_a+g*B_a, deg h=4, deg g=3, A_a,B_a affine")
    print("normalized affine poles alpha:", [float(x.midpoint()) for x in alphas])
    print("positive d parameters:", [float(x.midpoint()) for x in ds])
    print("norm-pair matrix principal 012 minor:", float(principal_minor.midpoint()))
    print("All divisions and nonzero tests used exact rational intervals.")


if __name__ == "__main__":
    main()
