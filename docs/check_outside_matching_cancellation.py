"""Exact outside-only matching cancellation; standard library arithmetic."""

from fractions import Fraction as Q
from itertools import combinations

from check_four_row_rational_quartic import (
    G, multiply, polynomial, monic_exact_division, trim,
)


def add(p, q):
    n = max(len(p), len(q))
    return trim([(p[j] if j < len(p) else G(0))
                 + (q[j] if j < len(q) else G(0)) for j in range(n)])


def scale(p, c):
    return trim([c * a for a in p])


def conjugate(p):
    return [a.conjugate() for a in p]


def norm(p):
    return multiply(p, conjugate(p))


def bracket(p, q):
    return trim([G(c.b) for c in multiply(p, conjugate(q))])


def main():
    rs = [Q(2), -Q(51, 22), Q(11, 3), -Q(1, 17)]
    zs = [G((1-r*r)/(1+r*r), 2*r/(1+r*r)) for r in rs]
    b = -sum(z.b for z in zs)
    c = -sum(z.a*z.b for z in zs)
    ws = [G((z.a*z.b+c)/(z.b+b), z.b+b) for z in zs]
    hs = [polynomial([ws[j]]+[zs[l] for l in range(4) if l != j])
          for j in range(4)]
    tau = [h[0].b for h in hs]
    norms = [norm(polynomial([z])) for z in zs]
    plus_one = [G(1), G(0), G(1)]
    assert add(scale(norms[0], 2117), scale(norms[1], -1851)) == scale(plus_one, 266)

    # Reanchor the original outside labels 0,1,2,3 at original row 4.
    rows = [conjugate(hs[3])]
    for j in range(3):
        shared = polynomial([zs[l] for l in range(4) if l not in (j, 3)])
        rows.append(monic_exact_division(multiply(hs[j], conjugate(hs[3])), norm(shared)))
    expected_t = [-tau[3]]+[tau[j]-tau[3] for j in range(3)]
    assert all(t != 0 for t in expected_t)
    for row, t in zip(rows, expected_t):
        assert len(row) == 5 and row[-1] == 1
        assert [a.b for a in row] == [t, Q(0), Q(0), Q(0), Q(0)]

    # Each matching acquires the same positive reanchoring factor.
    a1 = -tau[0]*(tau[1]-tau[2])
    a2 = -tau[1]*(tau[0]-tau[2])
    private_norm_square = norm(polynomial([ws[3]]))
    private_norm_square = multiply(private_norm_square, private_norm_square)
    m1 = multiply(bracket(rows[0], rows[1]), bracket(rows[2], rows[3]))
    m2 = multiply(bracket(rows[0], rows[2]), bracket(rows[1], rows[3]))
    assert m1 == scale(multiply(private_norm_square, multiply(norms[0], norms[3])), a1)
    assert m2 == scale(multiply(private_norm_square, multiply(norms[1], norms[3])), a2)
    lhs = add(scale(m1, Q(2117)/a1), scale(m2, -Q(1851)/a2))
    rhs = scale(multiply(private_norm_square, multiply(norms[3], plus_one)), 266)
    assert lhs == rhs

    # The new roots +/-i are absent from every old factor and its conjugate.
    roots = zs+ws
    assert all(z != G(0, sign) for z in roots for sign in (-1, 1))
    # Nor does any individual outside bracket contain either new root.
    for j, l in combinations(range(4), 2):
        d = bracket(rows[j], rows[l])
        for z in [G(0, 1), G(0, -1)]:
            value = G(0)
            for coef in reversed(d):
                value = value*z+coef
            assert value != 0
    print("Passed: four equal-degree constant-imaginary reanchored rows; "
          "two exact matching factorizations; new X^2+1 divisor of their "
          "fixed-coefficient sum; no new root in any individual bracket.")


if __name__ == "__main__":
    main()
