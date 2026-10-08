#!/usr/bin/env python3
"""Exact frame reconstruction fixtures, not a rational full K5 profile."""

from fractions import Fraction as F
from itertools import combinations

from check_four_row_rational_quartic import G, multiply, trim


def divide(a, b):
    a, b = G(a), G(b)
    norm = b.a*b.a+b.b*b.b
    assert norm
    c = a*b.conjugate()
    return G(c.a/norm, c.b/norm)


def evaluate(p, z):
    answer = G(0)
    for coefficient in reversed(p):
        answer = answer*z+coefficient
    return answer


def solve(matrix, rhs):
    """Gaussian elimination over Q(i), also returning the determinant."""
    n = len(matrix)
    a = [[G(x) for x in row]+[G(b)] for row, b in zip(matrix, rhs)]
    determinant = G(1)
    for j in range(n):
        pivot = next((k for k in range(j, n) if a[k][j] != 0), None)
        assert pivot is not None
        if pivot != j:
            a[j], a[pivot] = a[pivot], a[j]
            determinant = -determinant
        c = a[j][j]
        determinant = determinant*c
        a[j] = [divide(x, c) for x in a[j]]
        for k in range(n):
            if k != j:
                c = a[k][j]
                a[k] = [x-c*y for x, y in zip(a[k], a[j])]
    return [row[-1] for row in a], determinant


def bracket(h, g):
    return trim([G(c.b) for c in multiply([z.conjugate() for z in h], g)])


def contact_row(z, direction):
    """The scalar contact h(z)+direction*g(z)=0."""
    powers = [G(1)]
    for _ in range(4):
        powers.append(powers[-1]*z)
    return powers+[direction*x for x in powers[:4]]


def reconstruct(rows, indices):
    cols = [j for j in range(9) if j != 4]
    solution, determinant = solve(
        [[rows[k][j] for j in cols] for k in indices],
        [-rows[k][4] for k in indices])
    coefficients = solution[:4]+[G(1)]+solution[4:]
    return coefficients[:5], coefficients[5:], determinant


def residual(row, h, g):
    return sum((a*b for a, b in zip(row, h+g)), G(0))


def sylvester_frame(p, r):
    columns = []
    for j in range(4):
        columns.append([G(0)]*j+[-G(x) for x in r])
    for j in range(3):
        columns.append([G(0)]*j+[G(x) for x in p])
    matrix = [[column[k] if k < len(column) else G(0)
               for column in columns] for k in range(7)]
    solution, determinant = solve(matrix, [1]+[0]*6)
    assert all(c.b == 0 for c in solution)
    q, s = solution[:4], solution[4:]
    h = [G(a, b.a) for a, b in zip(p, q+[G(0)])]
    g = [G(a, b.a) for a, b in zip(r, s+[G(0)])]
    assert bracket(h, g) == [G(1)]
    return h, g, determinant


def main():
    # p*s-r*q=1 with q=r. These are genuine degree-(4,3) frames.
    p = [1, 0, 3, 0, 1]
    r = [0, 2, 0, 1]
    h, g, resultant = sylvester_frame(p, r)
    assert h == [G(1), G(0, 2), G(3), G(0, 1), G(1)]
    assert g == [G(0, 1), G(2), G(0, 1), G(1)]
    assert resultant in (G(1), G(-1))
    assert all(c.a.denominator == c.b.denominator == 1 for c in h+g)

    # Res(p,2*r)=2^4 Res(p,r). q is divided by two; s is unchanged.
    h2, g2, resultant2 = sylvester_frame(p, [2*a for a in r])
    assert resultant2 == 16*resultant
    assert any(c.b.denominator == 2 for c in h2+g2)
    assert all(abs(resultant2.a) % c.b.denominator == 0 for c in h2+g2)

    nodes = [G(k, F(k+1, 7)) for k in range(1, 11)]
    assert all(evaluate(g, z) != 0 for z in nodes)
    rows = [contact_row(z, -divide(evaluate(h, z), evaluate(g, z)))
            for z in nodes]
    # All 45 eight-contact subsets recover exactly the same rational frame.
    for indices in combinations(range(10), 8):
        hh, gg, determinant = reconstruct(rows, indices)
        assert determinant != 0
        assert hh == h and gg == g
        assert all(residual(row, hh, gg) == 0 for row in rows)

    # The two final contact residuals are indispensable.
    altered = [row[:] for row in rows]
    altered[8][5] = altered[8][5]+1
    hh, gg, _ = reconstruct(altered, range(8))
    assert residual(altered[8], hh, gg) != 0
    assert residual(altered[9], hh, gg) == 0

    # Ten contacts and full rank do not imply a constant frame determinant.
    bad_h = [G(0, 1), G(0), G(0), G(0), G(1)]
    bad_g = [G(0, 1), G(0), G(0), G(1)]
    bad_rows = [contact_row(z, -divide(evaluate(bad_h, z), evaluate(bad_g, z)))
                for z in nodes]
    hh, gg, _ = reconstruct(bad_rows, range(8))
    assert hh == bad_h and gg == bad_g
    assert all(residual(row, hh, gg) == 0 for row in bad_rows)
    assert bracket(hh, gg) == [G(0), G(0), G(0), G(-1), G(1)]

    print("PASS: every eight-contact subset reconstructs the same frame")
    print("PASS: remaining-contact and constant-determinant tests are independent")
    print("PASS: Sylvester reconstruction and nonunit resultant denominators")
    print("Scope: scalar interpolation fixtures only, not a rational K5 profile")


if __name__ == "__main__":
    main()
