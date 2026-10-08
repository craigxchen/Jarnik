#!/usr/bin/env python3
"""Bounded exact checks for the inert radical forced by cotangent cliques."""

from itertools import combinations
from math import gcd


def is_prime(n):
    return n > 1 and all(n % d for d in range(2, int(n ** 0.5) + 1))


def inert_primes(limit):
    return [p for p in range(3, limit + 1, 4) if is_prime(p)]


def compatible(a, b, L):
    d = a - b
    return d != 0 and (a * b + L * L) % d == 0


def maximal_cliques(vertices, adjacency):
    out = []

    def bronk(R, P, X):
        if not P and not X:
            out.append(tuple(sorted(R)))
            return
        # A pivot only reduces branching; deterministic choice is sufficient.
        u = next(iter(P | X), None)
        candidates = P - (adjacency[u] if u is not None else set())
        for v in list(candidates):
            bronk(R | {v}, P & adjacency[v], X & adjacency[v])
            P.remove(v)
            X.add(v)

    bronk(set(), set(vertices), set())
    return out


def required_exponent(p, m):
    e = 0
    power = 1
    while (p + 1) * power < m + 1:
        e += 1
        power *= p
    return e


def required_dyadic_exponent(m):
    e = 0
    power = 1
    while 4 * power < m + 1:
        e += 1
        power *= 2
    return e


def check_clique(A, L):
    m = len(A)
    assert all(compatible(a, b, L) for a, b in combinations(A, 2))
    for p in inert_primes(max(m, 3)):
        e = 0
        q = p
        while L % q == 0:
            e += 1
            q *= p
        assert e >= required_exponent(p, m)
    e2 = 0
    q = 2
    while L % q == 0:
        e2 += 1
        q *= 2
    assert e2 >= required_dyadic_exponent(m)


if __name__ == "__main__":
    all_seen = 0
    for L in range(1, 31):
        # This finite window is only a sanity check; divisor searches can be
        # substituted without changing the local theorem being tested.
        V = list(range(-60, 61))
        adj = {a: {b for b in V if b != a and compatible(a, b, L)} for a in V}
        cliques = maximal_cliques(V, adj)
        for A in cliques:
            check_clique(A, L)
        largest = max((len(A) for A in cliques), default=0)
        all_seen += len(cliques)
        print(f"L={L:2d}: maximal cliques={len(cliques):4d}, largest={largest}")
    assert all_seen > 0
    # The known sharp dyadic example is included explicitly.
    A = (-18, -6, -3, 0, 2, 6, 12)
    check_clique(A, 6)
    print("PASS: inert prime-power and dyadic exponent requirements")
