"""Exact simultaneous-alignment overlap divisibility and norm checks."""

from itertools import combinations
from math import factorial, gcd, prod
from random import Random

from check_endpoint_swap_content_and_overlap import (
    primitive_radius, reduce_denominator,
)
from check_least_radius_formula import (
    conj, exact_div, gcd_all, gcd_gaussian, lcm_all, mul, norm,
)


def product_gaussian(values):
    result = (1, 0)
    for value in values:
        result = mul(result, value)
    return result


def excess(B, row):
    A = reduce_denominator(conj(row))
    return exact_div(A, gcd_gaussian(B, A))


def check_case(x, rows, parameters):
    """Check the integer divisibility, not only its numerical upper bound."""
    T, B = x*x+1, (x, 1)
    assert x > 0 and x % 2 == 0
    assert len(set(parameters)) == len(parameters) > 0
    assert all(k > 0 and gcd(k, T) == 1 for k in parameters)
    assert all(d > 0 and v > 0 and gcd(d, v) == 1 for d, v in rows)
    assert len(set(rows)) == len(rows)
    r, K = len(rows), max(parameters)
    F = [gcd_gaussian(B, (d, 0)) for d, v in rows]
    coeffs = [(exact_div((d, 0), f), exact_div((x*v, v), f))
              for (d, v), f in zip(rows, F)]
    assert all(norm(gcd_gaussian(a, b)) == 1 for a, b in coeffs)
    Q = [excess(B, (d+x*v, v)) for d, v in rows]
    N = T*norm(lcm_all(Q))
    assert N == primitive_radius([(1, 0), B]+[(d+x*v, v) for d, v in rows])
    triples = list(combinations(range(r), 3))
    for triple in triples:
        subset_radius = T*norm(lcm_all([Q[i] for i in triple]))
        assert N % subset_radius == 0
    resultants = []
    for i, j in combinations(range(r), 2):
        di, vi = rows[i]
        dj, vj = rows[j]
        cross = di*vj-dj*vi
        assert cross != 0
        R = exact_div((x*cross, cross), mul(F[i], F[j]))
        resultants.append(R)
        assert norm(R)*norm(F[i])*norm(F[j]) == T*cross*cross
        assert norm(R) <= 4*norm(Q[i])*norm(Q[j])
    lhs = 1
    pair_products = [(1, 0)]*len(resultants)
    ordinary_contents = []
    for k in parameters:
        C = [(d+k*x*v, k*v) for d, v in rows]
        L = [exact_div(c, f) for c, f in zip(C, F)]
        Qp = [excess(B, c) for c in C]
        for (d, v), c, ell, q in zip(rows, C, L, Qp):
            gamma = gcd(d, k)
            ordinary_contents.append(gamma)
            assert gcd(*c) == gamma
            exact_div(ell, q)
            assert norm(q) >= 2*x*(k//gamma)*v
        overlap = prod(norm(q) for q in Qp)//norm(lcm_all(Qp))
        raw_overlap = prod(norm(ell) for ell in L)//norm(lcm_all(L))
        assert raw_overlap % overlap == 0
        Np = T*norm(lcm_all(Qp))
        assert Np == primitive_radius([(1, 0), B]+C)
        for triple in triples:
            subset_radius = T*norm(lcm_all([Qp[i] for i in triple]))
            assert Np % subset_radius == 0
        assert overlap*Np >= T*(2*x)**r*prod((k//gcd(d, k))*v for d, v in rows)
        lhs *= overlap
        for n, (i, j) in enumerate(combinations(range(r), 2)):
            pair_products[n] = mul(pair_products[n], gcd_gaussian(L[i], L[j]))
    fac = factorial(K-1)
    for R, pair_product in zip(resultants, pair_products):
        exact_div((fac*R[0], fac*R[1]), pair_product)
    rhs = fac**(r*(r-1))*prod(norm(R) for R in resultants)
    assert rhs % lhs == 0
    assert lhs <= (2*fac)**(r*(r-1))*prod(norm(q) for q in Q)**(r-1)
    assert lhs <= (2*fac)**(r*(r-1))*(N//T)**(r*(r-1))
    return max(ordinary_contents), lhs


def main():
    count = 0
    # Both nonconsecutive sets include nontrivial raw coordinate contents.
    for parameters in ([1, 2, 4, 8, 16], [2, 3, 7, 11]):
        check_case(18, [(2, 1), (3, 2), (125, 2), (650, 3)], parameters)
        count += 1
    # Exact common split-prime powers at source k=1 and target k=3.
    modulus = 13**3
    root = next(a for a in range(modulus) if (a*a+1) % modulus == 0)
    for target_k in (1, 3):
        x = 2
        first = target_k*(root-x) % modulus + modulus
        rows = [(first+j*modulus, 1) for j in range(4)]
        content, overlap = check_case(x, rows, [1, 2, 3, 4, 8, 16])
        assert overlap >= modulus**3
        target_excesses = [excess((x, 1), (d+target_k*x*v, target_k*v))
                           for d, v in rows]
        assert norm(gcd_all(target_excesses)) % modulus == 0
        target_overlap = (prod(norm(q) for q in target_excesses)
                          // norm(lcm_all(target_excesses)))
        assert target_overlap % (modulus**3) == 0
        count += 1
    rng = Random(713)
    for x in (2, 4, 6, 8, 10, 12, 18, 24):
        T = x*x+1
        for _ in range(100):
            rows, wanted = [], rng.randrange(3, 8)
            while len(rows) < wanted:
                d, v = rng.randrange(1, 3*T+1), rng.randrange(1, 31)
                if gcd(d, v) == 1 and (d, v) not in rows:
                    rows.append((d, v))
            parameters = [k for k in range(1, rng.randrange(3, 13))
                          if gcd(k, T) == 1 and rng.randrange(3) != 0]
            if not parameters:
                parameters = [1]
            check_case(x, rows, parameters)
            count += 1
    print('simultaneous alignment capacity cases:', count)


if __name__ == '__main__':
    main()
