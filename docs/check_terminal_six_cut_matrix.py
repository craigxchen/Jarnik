"""Exact 15-by-5 terminal coherent-cut matrix for six rows."""

from fractions import Fraction
from itertools import combinations, permutations
from math import gcd, lcm, prod


TRIANGLES = []
for first in combinations(range(6), 3):
    if 0 in first:
        TRIANGLES.append((first, tuple(i for i in range(6) if i not in first)))
BASIS_INDICES = (0, 1, 2, 4, 5)
BASIS = [TRIANGLES[i] for i in BASIS_INDICES]
CUTS = list(combinations(range(6), 2))


def add(a, b):
    n = max(len(a), len(b))
    out = [0] * n
    for i, x in enumerate(a):
        out[i] += x
    for i, x in enumerate(b):
        out[i] += x
    while len(out) > 1 and out[-1] == 0:
        out.pop()
    return tuple(out)


def mul(a, b):
    out = [0] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            out[i+j] += x*y
    while len(out) > 1 and out[-1] == 0:
        out.pop()
    return tuple(out)


def neg(a):
    return tuple(-x for x in a)


def b_polys():
    return [(0,), (1,), (2,), (3,), (4,), (0, 1)]


def bracket(i, j, inside, b):
    if i in inside and j in inside:
        return (0,)
    if i in inside:
        return (-1,)
    if j in inside:
        return (1,)
    return add(b[j], neg(b[i]))


def entry(partition, inside, b):
    value = (1,)
    for tri in partition:
        if len(set(tri) & inside) != 1:
            return (0,)
        a, c, d = tri
        for i, j in ((a, c), (c, d), (d, a)):
            value = mul(value, bracket(i, j, inside, b))
    return value


def matrix_polynomial():
    b = b_polys()
    return [[entry(part, set(cut), b) for part in BASIS] for cut in CUTS]


def evaluate_poly(poly, t):
    return sum(c*t**i for i, c in enumerate(poly))


def matrix_at(t):
    return [[evaluate_poly(a, t) for a in row] for row in matrix_polynomial()]


def rank(a):
    a = [[Fraction(x) for x in row] for row in a]
    r = 0
    for c in range(len(a[0]) if a else 0):
        p = next((i for i in range(r, len(a)) if a[i][c]), None)
        if p is None:
            continue
        a[r], a[p] = a[p], a[r]
        z = a[r][c]
        a[r] = [x/z for x in a[r]]
        for i in range(r+1, len(a)):
            if a[i][c]:
                z = a[i][c]
                a[i] = [x-z*y for x, y in zip(a[i], a[r])]
        r += 1
    return r


def primitive_kernel(a):
    a = [[Fraction(x) for x in row] for row in a]
    pivots, r = [], 0
    for c in range(len(a[0])):
        p = next((i for i in range(r, len(a)) if a[i][c]), None)
        if p is None:
            continue
        a[r], a[p] = a[p], a[r]
        z = a[r][c]
        a[r] = [x/z for x in a[r]]
        for i in range(len(a)):
            if i != r and a[i][c]:
                z = a[i][c]
                a[i] = [x-z*y for x, y in zip(a[i], a[r])]
        pivots.append(c)
        r += 1
    free = [j for j in range(len(a[0])) if j not in pivots]
    assert len(free) == 1
    v = [Fraction(0)] * len(a[0])
    v[free[0]] = 1
    for i, c in enumerate(pivots):
        v[c] = -a[i][free[0]]
    den = lcm(*(x.denominator for x in v))
    v = [int(x*den) for x in v]
    g = 0
    for x in v:
        g = gcd(g, abs(x))
    v = [x//g for x in v]
    if next(x for x in v if x) < 0:
        v = [-x for x in v]
    return v


def det4(a):
    out = 0
    for p in permutations(range(4)):
        inv = sum(p[i] > p[j] for i, j in combinations(range(4), 2))
        out += (-1)**inv * prod(a[i][p[i]] for i in range(4))
    return out


def maximal_minor_data(a):
    g, max_abs = 0, 0
    for rr in combinations(range(15), 4):
        for cc in combinations(range(5), 4):
            d = det4([[a[i][j] for j in cc] for i in rr])
            g = gcd(g, abs(d))
            max_abs = max(max_abs, abs(d))
    return g, max_abs


def primitive_formula(t):
    raw = (10*t-22, -6*t+18, 4*t-16, 2*t-6, 4-t)
    g = 0
    for x in raw:
        g = gcd(g, abs(x))
    out = tuple(x//g for x in raw)
    if next(x for x in out if x) < 0:
        out = tuple(-x for x in out)
    return out


def main():
    # In the gradient basis (G_A,G_B,G_C,G_D,G_E), these five triangle
    # products have the displayed unimodular coordinate matrix. The other
    # five complementary triangle products are integral combinations of them
    # by the exact identities in segre_gradient_arithmetic.md.
    triangle_to_gradient = [
        (0, 0, 0, 0, 1),
        (0, 0, 0, -1, 1),
        (0, 0, 1, -1, 0),
        (0, 1, 0, -1, 0),
        (-1, 1, 1, -1, 1),
    ]
    # Compute its determinant exactly to record saturation of this basis.
    det_basis = 0
    for p in permutations(range(5)):
        inv = sum(p[i] > p[j] for i, j in combinations(range(5), 2))
        det_basis += (-1)**inv * prod(triangle_to_gradient[i][p[i]]
                                      for i in range(5))
    assert abs(det_basis) == 1

    M = matrix_polynomial()
    q = ((-22, 10), (18, -6), (-16, 4), (-6, 2), (4, -1))
    # Exact polynomial identity: every row annihilates the kernel vector.
    for row in M:
        total = (0,)
        for a, b in zip(row, q):
            total = add(total, mul(a, b))
        assert total == (0,), total

    expected = {
        5: (1, 19600), 6: (4, 51984), 7: (1, 112896),
        8: (4, 215296), 9: (1, 374544), 10: (4, 608400),
        11: (1, 937024), 12: (4, 1382976),
    }
    for t, minor_data in expected.items():
        a = matrix_at(t)
        assert rank(a) == 4
        v = primitive_kernel(a)
        assert tuple(v) == primitive_formula(t), (t, v, primitive_formula(t))
        assert all(sum(x*y for x, y in zip(row, v)) == 0 for row in a)
        got = maximal_minor_data(a)
        assert got == minor_data, (t, got, minor_data)
        print(f"t={t}: kernel={v}, 4-minor-content={got[0]}, maxabs={got[1]}")
    print("PASS: unimodular triangle basis, symbolic 15x5 terminal matrix, kernel, and maximal-minor fixtures.")


if __name__ == "__main__":
    main()
