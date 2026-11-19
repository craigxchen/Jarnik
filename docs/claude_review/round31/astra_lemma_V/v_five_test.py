"""Certify the polynomial five-block model in Part V over Q(i)."""

from itertools import combinations

import sympy as sp


def norm(polynomial):
    return sp.expand(polynomial * sp.conjugate(polynomial))


def determinant(first, second):
    return sp.im(sp.expand(sp.conjugate(first) * second)).expand()


def fixed_nonzero_cofactor(first, second, modulus):
    cofactor = sp.cancel(determinant(first, second) / modulus)
    assert cofactor.is_Rational and cofactor != 0
    return cofactor


def verify_incidences(columns, edge, block, variable, expected):
    for mark in edge:
        remainder = sp.rem(norm(columns[mark]), block, variable)
        assert (remainder == 0) == expected


def main():
    x, parameter = sp.symbols("x T", real=True)
    marks = tuple(map(sp.Rational, ("3/2", "3", "1/3", "-32/15", "-5/54")))
    edges = tuple(combinations(range(len(marks)), 2))
    shared_roots = {(j, k): marks[j] * marks[k] for j, k in edges}
    assert len({abs(root) for root in shared_roots.values()}) == 10

    central = sp.expand(sp.prod(x + sp.I * shared_roots[0, j]
                                for j in range(1, len(marks))))
    auxiliary = ((4500480*x**2 + 1413120) / 1226641
                 + sp.I * (x**2 + sp.Rational(26758, 283071)*x + 1))
    assert determinant(central, auxiliary) == sp.expand((x - 1)**6)

    squared = [mark**2 for mark in marks]
    multiplier = -sp.Rational(418284581, 1009043070)
    columns = [
        sp.expand(central + multiplier * auxiliary * (u - squared[0]) * (squared[0] + x**2/u))
        for u in squared
    ]
    extra_root = sp.I * shared_roots[1, 2]
    assert all(sp.expand(columns[j].subs(x, extra_root)) == 0 for j in (1, 2))
    isotropic = {(0, 1), (0, 2), (0, 3), (0, 4), (1, 2)}
    for edge, signed_root in shared_roots.items():
        block = x**2 + signed_root**2
        j, k = edge
        cofactor = fixed_nonzero_cofactor(columns[j], columns[k], (x - 1)**6 * block)
        assert cofactor == multiplier * squared[0] * (1/squared[j] - 1/squared[k])
        verify_incidences(columns, edge, block, x, edge in isotropic)

    values_at_one = [column.subs(x, 1) for column in columns]
    scales = [sp.cancel(value / values_at_one[0]) for value in values_at_one]
    assert all(scale.is_Rational and scale != 0 for scale in scales)
    monic_columns = [sp.cancel(parameter**4 * column.subs(x, 1 + 1/parameter) / value).expand()
                     for column, value in zip(columns, values_at_one)]
    for column in monic_columns:
        polynomial = sp.Poly(column, parameter)
        assert polynomial.degree() == 4 and polynomial.LC() == 1
        assert sp.gcd(sp.re(column).expand(), sp.im(column).expand()) == 1

    for edge, signed_root in shared_roots.items():
        block = sp.expand(((parameter + 1)**2 + signed_root**2 * parameter**2)
                          / (1 + signed_root**2))
        root = sp.cancel(1 / (sp.I * signed_root - 1))
        assert sp.expand((parameter - root) * (parameter - sp.conjugate(root)) - block) == 0
        j, k = edge
        fixed_nonzero_cofactor(monic_columns[j], monic_columns[k], block)
        verify_incidences(monic_columns, edge, block, parameter, edge in isotropic)

    print("Ten signed determinants and twenty incidences verified in both parameter charts.")
    print("The isotropic blocks are 12, 13, 14, 15 and 23 (one-based labels).")
    print("The transformed columns are monic and primitive over Q[T].")


if __name__ == "__main__":
    main()
