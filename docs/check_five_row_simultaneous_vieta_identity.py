"""Exact symbolic certificate for the common-Q five-swap identity.

All polynomial calculations use rational coefficients.  Fixtures also
include an infinite new root and independent projective row rescalings.
No height bound or endpoint theorem is asserted.
"""

from fractions import Fraction
from itertools import combinations
from math import prod

from check_five_row_gradient_discriminant_core_count import (
    evaluate, independent_graphs, inverse,
)
from check_smaller_cut_relation_rank import (
    add, derivative, difference, mul, scale, variable,
)


ONE = {(0,) * 5: 1}


def affine(poly):
    out = {}
    for monomial, coefficient in poly.items():
        key = monomial[1::2]
        out[key] = out.get(key, 0) + coefficient
    return {key: value for key, value in out.items() if value}


def total(polys):
    out = {}
    for poly in polys:
        out = add(out, poly)
    return out


def polynomial_product(polys):
    out = ONE
    for poly in polys:
        out = mul(out, poly)
    return out


def cubic(i, xs):
    others = [x for j, x in enumerate(xs) if j != i]
    s1 = total(others)
    s2 = total(mul(a, b) for a, b in combinations(others, 2))
    s3 = total(polynomial_product(triple)
               for triple in combinations(others, 3))
    return total((scale(polynomial_product([xs[i]] * 3), 8),
                  scale(polynomial_product([s1, xs[i], xs[i]]), -6),
                  scale(mul(s2, xs[i]), 4),
                  scale(polynomial_product([s1] * 3), -2),
                  scale(mul(s1, s2), 8), scale(s3, -18)))


def symbolic_certificate(basis):
    xs = [variable(5, i) for i in range(5)]
    s1 = total(xs)
    fp = [polynomial_product(difference(5, i, j)
                             for j in range(5) if j != i)
          for i in range(5)]
    tests = 0
    for q in basis:
        lam = [derivative(q, i) for i in range(5)]
        a = [scale(derivative(lam[i], i), Fraction(1, 2))
             for i in range(5)]
        # Ward identities before imposing Q=0.
        assert not total(lam)
        assert total(mul(xs[i], lam[i]) for i in range(5)) == scale(q, 5)
        assert total(polynomial_product([xs[i], xs[i], lam[i]])
                     for i in range(5)) == scale(mul(s1, q), 2)
        for i in range(5):
            left = scale(mul(fp[i], lam[i]), 6)
            for j in range(5):
                linear = total((scale(xs[i], 2), scale(xs[j], 3),
                                scale(s1, -1)))
                left = add(left, scale(polynomial_product(
                    [fp[j], linear, a[j]]), -1))
            assert left == mul(cubic(i, xs), q)
            tests += 1
    return tests


def symbolic_diagonal_determinant(basis):
    matrix = [basis] + [
        [scale(derivative(derivative(q, i), i), Fraction(1, 2))
         for q in basis] for i in range(5)]
    states = {0: ONE}
    for row in matrix:
        next_states = {}
        for mask, partial in states.items():
            for j, entry in enumerate(row):
                if mask & (1 << j):
                    continue
                inversions = bin(mask >> (j + 1)).count("1")
                key = mask | (1 << j)
                term = scale(mul(partial, entry), (-1)**inversions)
                next_states[key] = add(next_states.get(key, {}), term)
        states = next_states
    vandermonde = polynomial_product(difference(5, j, i)
                                     for i, j in combinations(range(5), 2))
    assert states[63] == scale(mul(vandermonde, vandermonde), 6)


def bracket(v, w):
    return v[0] * w[1] - v[1] * w[0]


def cross_sum(old, new):
    return sum(Fraction(bracket(old[i], old[j]) * bracket(new[i], new[j]),
                        bracket(old[i], new[i]) * bracket(old[j], new[j]))
               for i, j in combinations(range(5), 2))


def check_zero(q, xs, basis):
    assert evaluate(q, xs) == 0
    lam = [evaluate(derivative(q, i), xs) for i in range(5)]
    if not all(lam):
        return False, False
    aa = [evaluate(derivative(derivative(q, i), i), xs) / 2
          for i in range(5)]
    assert all(isinstance(a, Fraction) for a in aa)
    old = [(Fraction(1), x) for x in xs]
    new = [(a, a * x - l) for a, x, l in zip(aa, xs, lam)]
    for i, (a, x, l) in enumerate(zip(aa, xs, lam)):
        # Homogeneous row quadratic evaluated at the new root, including a=0.
        linear = l - 2 * a * x
        constant = a * x * x - l * x
        X, Y = new[i]
        assert a * Y * Y + linear * X * Y + constant * X * X == 0
        assert bracket(old[i], new[i]) == -l
    assert cross_sum(old, new) == 6
    fp = [prod(xs[i] - xs[j] for j in range(5) if j != i)
          for i in range(5)]
    values = [l * f for l, f in zip(lam, fp)]
    A = (values[1] - values[0]) / (xs[1] - xs[0])
    B = values[0] - A * xs[0]
    assert all(v == A * x + B for v, x in zip(values, xs))
    u = [a / l for a, l in zip(aa, lam)]
    U = [sum(x**r * v for x, v in zip(xs, u)) for r in range(3)]
    assert A * U[1] + B * U[0] == 3 * A
    assert A * U[2] + B * U[1] == sum(xs) * A + 2 * B
    assert (U[1] - 3) * (U[1] - 2) == U[0] * (U[2] - sum(xs))
    # Reconstruct Q projectively from its five other roots.  The two
    # rows of the moment matrix cannot both vanish, since U_1 cannot
    # simultaneously equal 3 and 2.
    row = (U[1] - 3, U[0])
    if not any(row):
        row = (U[2] - sum(xs), U[1] - 2)
    new_A, new_B = row[1], -row[0]
    assert new_A or new_B
    target_a = [v * (new_A * x + new_B) / d
                for v, x, d in zip(u, xs, fp)]
    matrix = [[evaluate(p, xs) for p in basis]] + [
        [evaluate(derivative(derivative(p, i), i), xs) / 2
         for p in basis] for i in range(5)]
    inv = inverse(matrix)
    target = [Fraction(0)] + target_a
    coefficients = [sum(c * t for c, t in zip(row, target)) for row in inv]
    reconstructed = total(scale(p, c) for p, c in zip(basis, coefficients))
    factor = new_A / A if A else new_B / B
    assert reconstructed == scale(q, factor)
    # The rational expression is separately projective in every old/new row.
    old_scaled = [((i + 2) * v[0], (i + 2) * v[1])
                  for i, v in enumerate(old)]
    new_scaled = [((-1)**i * (i + 3) * v[0],
                   (-1)**i * (i + 3) * v[1]) for i, v in enumerate(new)]
    assert cross_sum(old_scaled, new_scaled) == 6
    # A common determinant-one projective transformation.
    transform = lambda v: (2 * v[0] + v[1], v[0] + v[1])
    assert cross_sum(list(map(transform, old)), list(map(transform, new))) == 6
    return True, any(a == 0 for a in aa)


def fixtures(basis):
    tests = infinity_tests = 0
    for source in ((0, 1, 2, 3, 4), (0, 1, 3, 7, 11),
                   (-3, 2, 5, 7, 12), (1, 2, 3, 5, 8),
                   (-4, -3, -2, -1, 1)):
        xs = tuple(map(Fraction, source))
        vals = [evaluate(q, xs) for q in basis]
        zeros = [add(scale(basis[0], vals[j]), scale(basis[j], -vals[0]))
                 for j in range(1, 6)]
        candidates = zeros[:]
        # Also impose a_i=0 to verify the homogeneous infinity case.
        for i in range(5):
            a = [evaluate(derivative(derivative(q, i), i), xs) for q in zeros]
            candidates.extend(add(scale(zeros[j], a[k]), scale(zeros[k], -a[j]))
                              for j, k in combinations(range(5), 2))
        for q in candidates:
            ok, infinity = check_zero(q, xs, basis)
            tests += ok
            infinity_tests += infinity
    assert tests >= 30 and infinity_tests >= 10
    return tests, infinity_tests


def main():
    basis = [affine(q) for _, q in independent_graphs(5, 2)]
    assert len(basis) == 6
    symbolic = symbolic_certificate(basis)
    symbolic_diagonal_determinant(basis)
    numeric, infinity = fixtures(basis)
    print(f"PASS: {symbolic} symbolic Hessian identities and exact diagonal determinant; "
          f"{numeric} exact five-swap/reconstruction fixtures, "
          f"{infinity} with an infinite root")


if __name__ == "__main__":
    main()
