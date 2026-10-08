"""Exact checks for actual cliques saturating the derivative depth bound."""
from itertools import combinations
from math import gcd, lcm, prod

from check_integer_cotangent_normalization import (
    edge_norm, edge_quotient, primitive_tuple,
)


def valuation(n, p):
    assert n
    v = 0
    while n % p == 0:
        n //= p
        v += 1
    return v


def root_minus_one(p, e):
    a = next(a for a in range(1, p) if (a*a+1) % p == 0)
    modulus = p
    for _ in range(1, e):
        a = next(a+c*modulus for c in range(p)
                 if ((a+c*modulus)**2+1) % (modulus*p) == 0)
        modulus *= p
    return a


def main():
    count = 0
    for k in range(4, 13):
        p = 13  # A split prime larger than every tested clique size.
        L = lcm(*range(1, k))
        for e in range(1, 5):
            P = p**e
            a0 = root_minus_one(p, e)
            a = next(a0+c*P for c in range(1, p+1)
                     if all(valuation((a0+(c+j)*P)**2+1, p) == e
                            for j in range(k)))
            assert P <= a <= (p+1)*P
            xs = [L*(a+j*P) for j in range(k)]
            edges = xs + [edge_quotient(x, y, L)
                          for x, y in combinations(xs, 2)]
            N = lcm(*(edge_norm(q, L) for q in edges))
            rows = primitive_tuple(xs, L)
            assert rows[0][0]**2 + rows[0][1]**2 == N
            assert valuation(N, p) == e
            for j, x in enumerate(xs):
                derivative = prod(x-y for t, y in enumerate(xs) if t != j)
                assert valuation(derivative, p) == (k-1)*e
                assert valuation(x*x+L*L, p) == e
                assert (x*x+L*L)**(k-1) % derivative == 0
            assert valuation(prod(x*x+L*L for x in xs), p) == k*e
            count += 1
    print(f"PASS: {count} actual positive cliques, all-edge lcm, Gaussian realization,")
    print("derivative saturation and resultant/conductor valuation mismatch.")


if __name__ == "__main__":
    main()
