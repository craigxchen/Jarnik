#!/usr/bin/env python3
"""Exact finite-cut square-quartic construction; no height conclusion."""

from fractions import Fraction
from itertools import combinations
from math import comb, prod

from check_gale_square_quartic_characterization import (
    interpolate,
    polynomial_multiply,
)


def valuation(value, prime):
    assert value
    exponent = 0
    while value % prime == 0:
        value //= prime
        exponent += 1
    return exponent


def split_primes(number, start=17):
    answer = []
    candidate = start
    while len(answer) < number:
        if candidate % 4 == 1 and all(
            candidate % divisor for divisor in range(2, int(candidate**0.5) + 1)
        ):
            answer.append(candidate)
        candidate += 1
    return answer


def lifted_root(prime, depth):
    root = next(r for r in range(prime) if (r * r + 1) % prime == 0)
    modulus = prime
    for _ in range(depth):
        root += modulus * (
            -((root * root + 1) // modulus) * pow(2 * root, -1, prime) % prime
        )
        modulus *= prime
    assert (root * root + 1) % modulus == 0
    return root


def add_congruence(value, modulus, residue, new_modulus):
    return value + modulus * ((residue - value) * pow(modulus, -1, new_modulus) % new_modulus)


def check_budgets():
    for n, expected in [(7, (154, 126, 35, 0, 63)), (8, (372, 406, 196, 35, 127))]:
        cuts = [subset for size in range(1, n // 2 + 1)
                for subset in combinations(range(n), size)
                if 2 * size != n or 0 in subset]
        terms = [sum(comb(len(subset), k) for subset in cuts) for k in range(1, 5)]
        assert tuple(terms + [terms[0] - terms[1] + terms[2] - terms[3]]) == expected
        assert len(cuts) == 2 ** (n - 1) - 1


def check_all_seven_cuts():
    n = 7
    cuts = [subset for size in range(1, 4) for subset in combinations(range(n), size)]
    primes = split_primes(len(cuts))
    nodes = [0] * n
    modulus = 1
    specifications = []
    for index, (subset, prime) in enumerate(zip(cuts, primes)):
        depth = 1 + (index % 3)
        root = lifted_root(prime, depth)
        next_modulus = prime ** (depth + 1)
        minority = {label: rank + 1 for rank, label in enumerate(subset)}
        ordinary = iter(r for r in range(prime) if (r * r + 1) % prime)
        for label in range(n):
            residue = root + prime**depth * minority[label] if label in minority else next(ordinary)
            nodes[label] = add_congruence(nodes[label], modulus, residue, next_modulus)
        modulus *= next_modulus
        specifications.append((subset, prime, depth))
    assert len(set(nodes)) == n
    norms = [1 + t * t for t in nodes]
    differences = {(i, j): nodes[j] - nodes[i] for i, j in combinations(range(n), 2)}
    checks = 0
    for subset, prime, depth in specifications:
        for i in range(n):
            assert valuation(norms[i], prime) == (depth if i in subset else 0)
            checks += 1
        for (i, j), difference in differences.items():
            assert valuation(difference, prime) == (depth if i in subset and j in subset else 0)
            checks += 1
        s = len(subset)
        r = sum(valuation(value, prime) for value in norms)
        d = sum(valuation(value, prime) for value in differences.values())
        t = sum(min(valuation(norms[i], prime) for i in triple) for triple in combinations(range(n), 3))
        assert (r, d, t) == (s * depth, comb(s, 2) * depth, comb(s, 3) * depth)
        assert r - d + t == depth
    derivatives = [prod(nodes[i] - nodes[j] for j in range(n) if i != j) for i in range(n)]
    weights = [Fraction(norms[i] ** 2, derivatives[i]) for i in range(n)]
    assert sum(weights) == 0
    assert sum(a * t for a, t in zip(weights, nodes)) == 0
    quartic = polynomial_multiply([1, 0, 1], [1, 0, 1])
    assert interpolate(nodes, [a * d for a, d in zip(weights, derivatives)]) == quartic
    return checks, max(t.bit_length() for t in nodes)


if __name__ == "__main__":
    check_budgets()
    checks, bits = check_all_seven_cuts()
    print("PASS: seven/eight resultant, discriminant, and higher-gcd cut budgets.")
    print(f"PASS: all 63 prescribed cuts, depths 1-3, {checks} exact valuations; nodes at most {bits} bits.")
    print("PASS: exact square-quartic interpolant and both seven-row Gale zero sums.")
