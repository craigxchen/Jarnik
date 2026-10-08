"""Exact checks for rational endpoint composition and orbit height."""

from itertools import combinations
from fractions import Fraction
from math import gcd

from check_endpoint_swap_content_and_overlap import all_edge_radius, primitive_radius


def rows(a, b, interiors, p=1, q=1):
    return [(1, 0), (a, b)] + [
        (q*d+p*a*v, p*b*v) for d, v in interiors
    ]


def radius(a, b, interiors, p=1, q=1):
    values = rows(a, b, interiors, p, q)
    n = all_edge_radius(values)
    assert n == primitive_radius(values)
    return n


def matmul(A, B):
    return tuple(tuple(sum(A[i][k]*B[k][j] for k in range(2))
                       for j in range(2)) for i in range(2))


def fixing_matrix(a, b, p, q):
    return ((q, Fraction(a*(p-q), b)), (0, p))


def swapping_matrix(a, b, p, q):
    T = a*a+b*b
    return ((q*a*b, p*T-q*a*a), (q*b*b, -q*a*b))


def check_orbit():
    fixtures = [
        (1, 1, [(1, 1), (2, 1), (6, 1)]),
        (2, 1, [(1, 1), (3, 2), (5, 1)]),
        (5, 2, [(1, 1), (3, 2), (7, 1)]),
        (7, 3, [(2, 1), (5, 2), (9, 1)]),
    ]
    params = [(1, 1), (2, 1), (1, 2), (3, 2), (2, 3), (5, 4)]
    relative_checks = 0
    reflection_checks = 0
    for a, b, interiors in fixtures:
        assert gcd(a, b) == 1
        assert all(gcd(d, v) == 1 for d, v in interiors)
        assert len({(d, v) for d, v in interiors}) == len(interiors)
        original = radius(a, b, interiors)
        T = a*a+b*b
        values = {}
        for p, q in params:
            n = radius(a, b, interiors, p, q)
            values[p, q] = n
            for d, v in interiors:
                # Exact scalar content bound, including primes dividing b.
                content = gcd(q*d+p*a*v, p*b*v)
                gamma = gcd(p, d)*gcd(q, v)
                assert gamma == gcd(q*d, p*v)
                assert content <= b*gamma
                assert gcd(d+a*v, b*v) <= b

            # Independent matrix multiplication confirms the group law
            # and the swap-to-fixing conversion with exact fractions.
            r, s = 3, 2
            assert matmul(fixing_matrix(a, b, r, s),
                          fixing_matrix(a, b, p, q)) == fixing_matrix(a, b, p*r, q*s)
            product = matmul(swapping_matrix(a, b, 1, 1),
                             swapping_matrix(a, b, p, q))
            expected = fixing_matrix(a, b, p, q)
            assert product == tuple(tuple(b*b*T*entry for entry in row)
                                    for row in expected)

        for (p, q), (r, s) in combinations(params, 2):
            g = gcd(r*q, s*p)
            P, Q = r*q//g, s*p//g
            assert a*a*(P+Q)**2 <= 4*b**4*values[p, q]*values[r, s]
            relative_checks += 1

        for (di, vi), (dj, vj) in combinations(interiors, 2):
            p, q = di*dj, T*vi*vj
            g = gcd(p, q)
            p, q = p//g, q//g
            pair = [(di, vi), (dj, vj)]
            assert radius(a, b, pair, p, q) == radius(a, b, pair)
            # G_(p/q)(H_i) and the conformal reflection of H_j
            # represent the same rational phase, and conversely.
            ci = (q*di+p*a*vi, p*b*vi)
            cj = (q*dj+p*a*vj, p*b*vj)
            ref_i = (T*vi+a*di, b*di)
            ref_j = (T*vj+a*dj, b*dj)
            assert ci[0]*ref_j[1] == ci[1]*ref_j[0]
            assert cj[0]*ref_i[1] == cj[1]*ref_i[0]
            reflection_checks += 1

        assert original == values[1, 1]

    # The edge to 2+i bounds p+q by 15 whenever the full radius is <=125.
    # Exhaust the resulting finite set, using the full all-edge lcm.
    interiors = [(1, 1), (2, 1), (6, 1)]
    assert radius(1, 1, interiors) == 125
    checked = 0
    minima = []
    for p in range(1, 15):
        for q in range(1, 16-p):
            if gcd(p, q) != 1:
                continue
            checked += 1
            n = radius(1, 1, interiors, p, q)
            if n <= 125:
                minima.append((n, p, q))
    assert checked == 71
    assert minima == [(125, 1, 1)]
    assert radius(1, 1, interiors, 3, 1) == 425
    return relative_checks, reflection_checks, checked


if __name__ == "__main__":
    print("relative height, pair reciprocity, no-descent:", check_orbit())
