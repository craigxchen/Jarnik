"""Exact checks for the polynomial families in Part VII. Uses SymPy over QQ(i)."""

from itertools import combinations
from functools import reduce
import sympy as sp


t = sp.Symbol("t", real=True)
field = sp.QQ_I


def poly(expression):
    return sp.Poly(expression, t, domain=field)


def bar(f):
    return poly(sp.conjugate(f.as_expr()))


def dimensions(marks):
    one = poly(1)
    parent = marks[0] * bar(marks[0])
    assert all(z * bar(z) == parent for z in marks)
    assert reduce(sp.gcd, marks).monic() == one
    degrees = {z.degree() for z in marks}
    assert len(degrees) == 1
    degree = degrees.pop()
    distances = [(z - w).degree() for z, w in combinations(marks, 2)]
    assert all(k >= 0 for k in distances)
    return degree, degree - max(distances)


def pell_blocks(n, slopes):
    blocks = []
    for slope in slopes:
        x, y = poly(1), poly(slope)
        for _ in range(n):
            x, y = poly(t*t+1)*x + poly(t)*y, poly(t)*x + y
        blocks.append(x + poly(sp.I)*y)
    return blocks


def words(blocks, edges):
    beta = poly(t + 2*sp.I)
    first = beta * sp.prod(bar(g) for g in blocks)
    return [first] + [
        bar(beta) * sp.prod(g if j in edge else bar(g)
                            for j, g in enumerate(blocks))
        for edge in edges
    ]


def main():
    triangle = [(0, 1), (0, 2), (1, 2)]
    for n in [1, 2, 3, 5]:
        for r, edges in [(3, triangle), (4, triangle + [(0, 3)])]:
            data = words(pell_blocks(n, [1, 3, 5, 8][:r]), edges)
            degree, depth = dimensions(data)
            assert (degree, depth) == (2*r*n+1, 4*n)
            print(f"marks={r+1} n={n} D={degree} h={depth}")
    for d in [1, 3, 6]:
        blocks = [poly(t**d-a+sp.I) for a in [1, 3, 8]]
        data = [sp.prod(g if k == j else bar(g) for k, g in enumerate(blocks))
                for j in range(3)]
        assert dimensions(data) == (3*d, 2*d)
        print(f"marks=3 d={d} D={3*d} h={2*d}")


if __name__ == "__main__":
    main()
