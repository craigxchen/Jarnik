#!/usr/bin/env python3
"""Exact certificate for one member of the vertical triangle family."""

from fractions import Fraction as F


def add(z, w):
    return (z[0] + w[0], z[1] + w[1])


def mul(z, w):
    return (z[0] * w[0] - z[1] * w[1], z[0] * w[1] + z[1] * w[0])


def product(values):
    out = (F(1), F(0))
    for value in values:
        out = mul(out, value)
    return out


def elementary_symmetric(roots, degree):
    # Small exact subset expansion, sufficient for degree four.
    from itertools import combinations

    total = (F(0), F(0))
    for indices in combinations(range(len(roots)), degree):
        total = add(total, product(roots[i] for i in indices))
    return total


def poly_mul(a, b):
    out = [F(0)] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            out[i + j] += x * y
    return out


def norm_linear(root):
    # Norm(T-root) = (T-Re(root))^2 + Im(root)^2, ascending coefficients.
    a, b = root
    return [a * a + b * b, -2 * a, F(1)]


def main():
    # Stereographic parameters u=1/2, v=1/4.
    u, v = F(1, 2), F(1, 4)
    den = 1 - u * u - v * v
    s = (1 + u * u + v * v) / den
    x = 2 * u / den
    y = 2 * v / den
    assert s * s - x * x - y * y == 1

    p, q, r = s * s, -x * x, -y * y
    X = s * x * y
    assert p + q + r == 1
    assert p * q * r == X * X

    roots = [
        (F(0), F(1)),
        (X, p - 1),
        (X, q - 1),
        (X, r - 1),
        (X * (1 + 1 / p), p),
        (X * (1 + 1 / q), q),
        (X * (1 + 1 / r), r),
    ]
    assert len(set(roots)) == 7
    for i, (a, b) in enumerate(roots):
        assert b != 0
        for c, d in roots[i + 1 :]:
            assert not (a == c and b == -d)

    rows = []
    # The row opposite t uses the two edge roots indexed by the other
    # parameters and has the private root indexed by t.
    row_data = {}
    for first, second, missing in ((p, q, r), (p, r, q), (q, r, p)):
        row = [
            (F(0), F(1)),
            (X, first - 1),
            (X, second - 1),
            (X * (1 + 1 / missing), missing),
        ]
        e1, e2, e3, e4 = (elementary_symmetric(row, k) for k in range(1, 5))
        assert e1[1] == e2[1] == e3[1] == 0
        expected_e1 = X * (3 + 1 / missing)
        expected_e2 = 1 + missing * missing + first * second + 3 * X * X
        expected_e3 = X * (
            1
            + missing * missing
            + (missing + 1) / missing
            + first * second * (missing * missing - 1) / missing
        )
        expected_real_e4 = X * X * (3 + 1 / missing) + missing * missing
        assert e1[0] == expected_e1
        assert e2[0] == expected_e2
        assert e3[0] == expected_e3
        expected_imaginary_constant = (
            X
            * (missing * missing - 1)
            * (X * X + missing * missing)
            / (missing * missing)
        )
        assert e4[1] == expected_imaginary_constant != 0
        assert e4[0] == expected_real_e4
        row_data[missing] = {
            "e1": e1[0],
            "e2": e2[0],
            "e3": e3[0],
            "real_e4": e4[0],
            "imag_e4": e4[1],
        }
        rows.append((e1, e2, e3, e4))

    # For rows missing a,b, their shared roots are i and the edge root
    # indexed by the third parameter c. The determinant is exactly
    # (Y_b-Y_a) Norm(T-i) Norm(T-(X+i(c-1))).
    for missing_a, missing_b, shared in ((p, q, r), (p, r, q), (q, r, p)):
        arow, brow = row_data[missing_a], row_data[missing_b]
        ya, yb = arow["imag_e4"], brow["imag_e4"]
        real_a = [arow["real_e4"], -arow["e3"], arow["e2"], -arow["e1"], F(1)]
        real_b = [brow["real_e4"], -brow["e3"], brow["e2"], -brow["e1"], F(1)]
        determinant = [yb * ca - ya * cb for ca, cb in zip(real_a, real_b)]
        shared_norm = poly_mul(
            norm_linear((F(0), F(1))), norm_linear((X, shared - 1))
        )
        expected = [(yb - ya) * coeff for coeff in shared_norm]
        assert determinant == expected

    assert [(row[3][1]) for row in rows] == [
        F(-53430249600, 2357947691),
        F(31636332000, 2357947691),
        F(1605496422400, 49516901511),
    ]
    print("u,v:", u, v)
    print("s,x,y:", s, x, y)
    print("p,q,r,X:", p, q, r, X)
    print("seven roots:")
    for root in roots:
        print(" ", root)
    print("row constant imaginary parts:", *(row[3][1] for row in rows))
    print("exact checks passed")


if __name__ == "__main__":
    main()
