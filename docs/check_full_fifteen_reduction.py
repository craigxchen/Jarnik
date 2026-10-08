"""Exact checks of the full-fifteen parameterization, not of its solvability.

Uses rational inputs and the standard-library Gaussian arithmetic from the
earlier quartic certificate. No input here is claimed to pass the six
remaining divisibility conditions.
"""

from fractions import Fraction as Q
from itertools import combinations
from random import Random

from check_four_row_rational_quartic import G, multiply, polynomial, trim


def inv(z):
    z = G(z)
    norm = z.a * z.a + z.b * z.b
    return G(z.a / norm, -z.b / norm)


def add(p, q):
    return trim([(p[j] if j < len(p) else G())
                 + (q[j] if j < len(q) else G())
                 for j in range(max(len(p), len(q)))])


def scale(p, c):
    return trim([a * c for a in p])


def sub(p, q):
    return add(p, scale(q, -1))


def conjugate(p):
    return [c.conjugate() for c in p]


def divmod_poly(p, q):
    p, q = trim(p), trim(q)
    assert q != [G()]
    result = [G()] * max(1, len(p) - len(q) + 1)
    while p != [G()] and len(p) >= len(q):
        shift, c = len(p) - len(q), p[-1] * inv(q[-1])
        result[shift] = c
        p = sub(p, [G()] * shift + scale(q, c))
    return trim(result), p


def exact_divide(p, q):
    result, remainder = divmod_poly(p, q)
    assert remainder == [G()]
    return result


def evaluate(p, z):
    result = G()
    for a in reversed(p):
        result = result * z + a
    return result


def product(polys):
    result = [G(1)]
    for p in polys:
        result = multiply(result, p)
    return result


def construct(roots, constants):
    eye, t = G(0, 1), [G(), G(1)]
    nodes = roots + [z.conjugate() for z in roots]
    values = [-inv(z-eye) for z in roots]
    values += [v.conjugate() for v in values]
    K = [G()]
    for j, (z, v) in enumerate(zip(nodes, values)):
        basis = polynomial(nodes[:j] + nodes[j+1:])
        K = add(K, scale(basis, v * inv(evaluate(basis, z))))
    assert len(K) <= 8 and all(a.b == 0 for a in K)
    c = K[7] if len(K) == 8 else G()
    quadratics = [polynomial([z, z.conjugate()]) for z in roots]
    P = [[G(d), c] for d in constants]
    Z = polynomial(roots)
    M = exact_divide(add(multiply([-eye, G(1)], K), [G(1)]), Z)
    assert len(M) <= 5
    bezout = sub(multiply([eye, G(1)], multiply(Z, M)),
                  multiply([-eye, G(1)], conjugate(multiply(Z, M))))
    assert bezout == [2 * eye]
    residuals, rows = [], []
    for j, z in enumerate(roots):
        others = [quadratics[k] for k in range(4) if k != j]
        A = sub(K, multiply(P[j], product(others)))
        assert len(A) <= 7
        F = add(multiply([G(1), G(), G(1)], A), t)
        W = sub(multiply([-z, G(1)], M),
                multiply(multiply([-eye, G(1)], P[j]),
                         polynomial([roots[k].conjugate()
                                     for k in range(4) if k != j])))
        assert len(W) <= 5
        reconstructed = multiply([eye, G(1)],
                                 multiply(polynomial([roots[k] for k in range(4)
                                                      if k != j]), W))
        assert reconstructed == add(F, [eye])
        assert evaluate(F, -eye) == -eye
        assert all(evaluate(F, roots[k]) == -eye for k in range(4) if k != j)
        rows.append(F)
        residuals.append(W)
    for j, k in combinations(range(4), 2):
        E = sub(multiply(P[k], quadratics[j]), multiply(P[j], quadratics[k]))
        assert len(E) <= 3
        predicted = multiply([G(1), G(), G(1)],
                             multiply(E, product([quadratics[l] for l in range(4)
                                                  if l not in (j, k)])))
        assert sub(rows[j], rows[k]) == predicted
    return K, rows


def main():
    rng = Random(20260910)
    for _ in range(24):
        roots = []
        while len(roots) < 4:
            z = G(Q(rng.randrange(-11, 12), rng.randrange(1, 5)),
                  Q(rng.randrange(1, 10), rng.randrange(1, 5)))
            if z in (G(0, 1), G(0, -1)) or z in roots:
                continue
            roots.append(z)
        construct(roots, [Q(rng.randrange(-5, 6), 3) for _ in range(4)])
    for a, b, d, e, p, q in [(2, 3, 4, -2, 5, 7), (1, 2, 3, 4, -2, 1)]:
        K, rows = construct([G(a,b), G(-a,b), G(d,e), G(-d,e)], [p,-p,q,-q])
        assert all(K[j] == 0 for j in range(0, len(K), 2))
        for j, k in [(0, 1), (2, 3)]:
            assert rows[k] == [a * (-1 if n % 2 == 0 else 1)
                               for n, a in enumerate(rows[j])]
    print("Passed: 24 rational inputs and two reflection slices; real CRT, Bezout, degree cancellations, four row factorizations, and six exact pair differences. Solvability is not asserted.")


if __name__ == '__main__':
    main()
