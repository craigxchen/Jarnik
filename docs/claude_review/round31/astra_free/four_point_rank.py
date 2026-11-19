"""Verify dominance of the rational four-point family by a modular Jacobian."""
from itertools import combinations

import sympy as sp


def run():
    modulus = 101
    imaginary_unit = 10
    assert imaginary_unit**2 % modulus == modulus - 1
    x, m, s = variables = sp.symbols("x m s")
    at_point = {x: -3, m: 3, s: 7}
    y = x * (1 - x) / (m**2 + x)
    z = 1 - x - y
    u = m * y

    def reduce(expression):
        value = sp.Rational(expression.subs(at_point))
        return int(value.p) * pow(int(value.q), -1, modulus) % modulus

    def make_quartic(a, b, c):
        return [(0, 1), (u, a - 1), (u, b - 1), (u * (1 + c) / c, c)]

    matrix = []
    for a, b, c in ((x, y, z), (x, z, y), (y, z, x)):
        roots = make_quartic(a, b, c)
        row = [0, 0, 0]
        for re, im in roots:
            for sign in (1, -1):
                factor = s - re - sign * imaginary_unit * im
                inverse = pow(reduce(factor), -1, modulus)
                for k, var in enumerate(variables):
                    row[k] += sign * reduce(sp.diff(factor, var)) * inverse
        matrix.append([v % modulus for v in row])
    determinant = int(sp.Matrix(matrix).det()) % modulus
    assert determinant
    print(f"Logarithmic Jacobian at x=-3, m=3, s=7 modulo {modulus}: {matrix}")
    print(f"Determinant modulo {modulus}: {determinant}")

    all_roots = set()
    all_roots.add((sp.Integer(0), sp.Integer(1)))
    for v in (x, y, z):
        all_roots.add((u, v - 1))
        all_roots.add((u * (1 + v) / v, v))
    values = [(sp.sympify(re).subs(at_point), sp.sympify(im).subs(at_point)) for re, im in all_roots]
    assert all(im for re, im in values)
    assert all((a, abs(b)) != (c, abs(d)) for (a, b), (c, d) in combinations(values, 2))
    print(f"Seven roots at this specialization: {sorted(values)}")


if __name__ == "__main__":
    run()
