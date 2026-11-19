"""Exact min-cost certificates for the five-coordinate grid examples in Part II."""

from fractions import Fraction
from itertools import product
from math import comb, prod

import numpy as np
from scipy.optimize import linear_sum_assignment


def total_degree(k, d):
    if k == 0:
        yield ()
        return
    for a in range(d + 1):
        for tail in total_degree(k - 1, d - a):
            yield (a,) + tail


def minimum(cost):
    rows, cols = linear_sum_assignment(cost)
    assert np.array_equal(rows, np.arange(len(rows)))
    assigned = cost[rows, cols]
    # Checking integer dual feasibility independently proves optimality of the
    # permutation found by scipy. Directed slack along i->j is the change in
    # row i's cost when it takes the column assigned to row j.
    edges = cost[:, cols] - assigned[:, None]
    potential = np.zeros(len(rows), dtype=np.int64)
    for _ in range(len(rows) - 1):
        relaxed = np.minimum(potential, (potential[:, None] + edges).min(axis=0))
        if np.array_equal(relaxed, potential):
            break
        potential = relaxed
    assert np.all(potential[None, :] <= potential[:, None] + edges)
    return int(assigned.sum())


def raw_exponent(rows, exponents, bits):
    valuations = np.array([[bits[j] for j in row] for row in rows], dtype=np.int64)
    exponents = np.array(exponents, dtype=np.int64)
    cost = valuations @ exponents.T
    degree = int(exponents.sum())
    return degree + minimum(cost) + minimum(-cost)


def shell_design(d):
    k = 5
    exponents = list(total_degree(k, d))
    new = 5 * (d + 1)
    node = lambda slot, level: 5 * level + slot
    special = [(d, 0, 0, 0, 0), (d - 1, 1, 0, 0, 0)]
    rows = []
    factors = []
    for alpha in exponents:
        if alpha in special:
            continue
        row = tuple(node(h, x) for h, x in enumerate(alpha))
        rows.append(row)
        for h, x in enumerate(alpha):
            for t in range(x):
                factors.append((row[h], node(h, t)))
    tail = tuple(node(h, 0) for h in (2, 3, 4))
    rows.extend([(new, new + 1) + tail, (new + 1, new) + tail])
    for j in (new, new + 1):
        for t in range(d - 1):
            factors.append((j, node(0, t)))
    factors.append((new, new + 1))
    factors.append((new, new + 1, node(0, d - 1), node(1, 0)))
    assert len(factors) == sum(map(sum, exponents))
    return exponents, rows, factors


def whole_polynomial_exponent(factors, bits):
    return sum(min(bits[j] for j in f) + 1 - max(bits[j] for j in f)
               for f in factors)


def walsh_cuts():
    return [[(x & v).bit_count() % 2 for x in range(32)] for v in range(1, 32)]


def exact_cyclic_determinant():
    from sympy import I, Matrix, Symbol, expand, resultant
    # Rational unit-circle parameters j/24 have total angular width < 10/24.
    t = 24
    radius = prod(t*t+j*j for j in range(6))
    z = [expand((t+I*j)**2 * (radius//(t*t+j*j))) for j in range(6)]
    assert all(expand(x*x.conjugate()) == radius**2 for x in z)
    matrix = Matrix([[1] + [z[(i+h)%6] for h in range(5)] for i in range(6)])
    determinant = expand(matrix.det(method='domain-ge'))
    x = Symbol('x')
    factors = resultant(sum(x**j for j in range(6)), sum(z[j]*x**j for j in range(6)), x)
    assert determinant != 0 and expand(determinant + factors) == 0
    print("Six integer circle points: Fourier product agrees up to sign with the nonzero determinant.")


def main():
    exact_cyclic_determinant()
    exponents = list(total_degree(5, 1))
    marks = [0, 1, 2, 4, 8, 16]
    cycle = [tuple(marks[(i + h) % 6] for h in range(5)) for i in range(6)]
    for pattern in product(range(2), repeat=6):
        bits = [0] * 32
        for j, b in zip(marks, pattern):
            bits[j] = b
        assert raw_exponent(cycle, exponents, bits) == 5 * int(len(set(pattern)) == 1)
    assert sum(raw_exponent(cycle, exponents, bits) for bits in walsh_cuts()) == 0
    print('Cyclic degree one: all 64 binary label patterns certified.')
    print('Raw cyclic kappa on 31 nonconstant Walsh columns: 0.')
    print('d    n    V    raw kappa     polynomial kappa')
    for d in (1, 2, 3, 4):
        exponents, rows, factors = shell_design(d)
        V = len(factors)
        assert V == 5 * comb(5 + d, 6)
        raw = after = 0
        for bits in walsh_cuts():
            local = raw_exponent(rows, exponents, bits)
            factored = whole_polynomial_exponent(factors, bits)
            assert local <= factored
            raw += local
            after += factored
        print(d, len(rows), V, Fraction(raw, 31*V), Fraction(after, 31*V))


if __name__ == '__main__':
    main()
