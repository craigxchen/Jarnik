#!/usr/bin/env python3
"""Exact checks of the even-minor formula and Dodgson identity for k<=4."""

from itertools import combinations
from fractions import Fraction
from math import isqrt


def det(a: list[list[int]]) -> int:
    n = len(a)
    if n == 0:
        return 1
    b = [row[:] for row in a]
    sign, prev = 1, 1
    for j in range(n - 1):
        if b[j][j] == 0:
            q = next((q for q in range(j + 1, n) if b[q][j]), None)
            if q is None:
                return 0
            b[j], b[q] = b[q], b[j]
            sign *= -1
        pivot = b[j][j]
        for r in range(j + 1, n):
            for c in range(j + 1, n):
                b[r][c] = (b[r][c] * pivot - b[r][j] * b[j][c]) // prev
        prev = pivot
    return sign * b[-1][-1]


def mul(z: tuple[int, int], w: tuple[int, int]) -> tuple[int, int]:
    return z[0] * w[0] - z[1] * w[1], z[0] * w[1] + z[1] * w[0]


def conj(z: tuple[int, int]) -> tuple[int, int]:
    return z[0], -z[1]


def scale(a: int, z: tuple[int, int]) -> tuple[int, int]:
    return a * z[0], a * z[1]


def norm(z: tuple[int, int]) -> int:
    return z[0] * z[0] + z[1] * z[1]


def basis(k: int):
    out = [(0, 0)]
    for j in range(1, k + 1):
        out += [(j, 0), (j - 1, 1)]
    return out


def matrix(points, columns):
    return [[x**a * y**b for a, b in columns] for x, y in points]


def omit(a, row, col):
    return [[v for j, v in enumerate(r) if j != col]
            for i, r in enumerate(a) if i != row]


def circle_points(n: int):
    pts = []
    for x in range(-isqrt(n), isqrt(n) + 1):
        y2 = n - x * x
        y = isqrt(max(0, y2))
        if y * y == y2:
            pts.extend([(x, y), (x, -y)] if y else [(x, 0)])
    return sorted(set(pts), key=lambda z: __import__("math").atan2(z[1], z[0]))


def main():
    n = 1105
    pts = circle_points(n)
    assert len(pts) >= 9
    units = ((1, 0), (-1, 0), (0, 1), (0, -1))
    for k in range(1, 5):
        seen = set()
        for start in range(min(5, len(pts) - 2 * k)):
            rows = pts[start:start + 2 * k]
            low = basis(k - 1)
            x = det(matrix(rows, low + [(k, 0)]))
            y = det(matrix(rows, low + [(k - 1, 1)]))
            e = (x, y)
            v = (1, 0)
            p = (1, 0)
            for z in rows:
                p = mul(p, z)
            for i, j in combinations(range(2 * k), 2):
                v = mul(v, (rows[j][0] - rows[i][0], rows[j][1] - rows[i][1]))
            assert norm(e) * 2 ** (2 * k * (k - 1)) * n ** (k * (k - 1)) == norm(v)
            exact_left = scale(2 ** (k * (k - 1)), e)
            for _ in range(k - 1):
                exact_left = mul(exact_left, p)
            epsilon = ((1, 0), (0, 1), (-1, 0), (0, -1))[(k - 1) % 4]
            exact_right = mul(epsilon, scale(n ** (k * (k - 1) // 2), v))
            assert exact_left == exact_right
            if k == 2:
                def difference(j, i):
                    return rows[j][0] - rows[i][0], rows[j][1] - rows[i][1]
                pairing_u = mul(difference(1, 0), difference(3, 2))
                pairing_v = mul(difference(2, 0), difference(3, 1))
                ratio_numerator = mul(pairing_v, conj(pairing_u))
                assert ratio_numerator[1] == 0
                ratio = Fraction(ratio_numerator[0], norm(pairing_u))
                alpha, beta = ratio.numerator, ratio.denominator
                assert alpha > beta > 0
                assert all(t % beta == 0 for t in pairing_u)
                gamma = tuple(t // beta for t in pairing_u)
                s_numerator = norm(gamma) * alpha * beta * (alpha - beta)
                assert s_numerator % (4 * n) == 0
                s = s_numerator // (4 * n)
                assert e == mul((0, 1), scale(s, gamma))
            left = scale(n**k, e)
            candidates = [mul(u, mul(conj(e), p)) for u in units]
            assert left in candidates
            seen.add(units[candidates.index(left)])
        assert seen == {(-1, 0)}

        rows = pts[:2 * k + 1]
        a = matrix(rows, basis(k))
        full = det(a)
        inner = det([row[:-2] for row in a[1:-1]])
        lc = det(omit(a, 0, 2 * k - 1))
        ld = det(omit(a, 0, 2 * k))
        rc = det(omit(a, 2 * k, 2 * k - 1))
        rd = det(omit(a, 2 * k, 2 * k))
        assert full * inner == lc * rd - ld * rc
    print("checked k=1..4 exact even-minor formulas, norms, phases, and Dodgson identities")


if __name__ == "__main__":
    main()
