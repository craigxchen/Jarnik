#!/usr/bin/env python3
"""Independent support/rank audit of the phase-coordinate Newton formulas."""

from itertools import combinations, permutations
from math import comb


N = 8
TRIPLES = list(combinations(range(N), 3))
CUTS = [s for k in range(1, 5) for s in combinations(range(N), k)
        if k != 4 or 0 in s]


def sign(permutation):
    return (-1) ** sum(permutation[i] > permutation[j]
                       for i in range(len(permutation))
                       for j in range(i + 1, len(permutation)))


def coordinate(triple, raw):
    complement = [i for i in range(N) if i not in triple]
    answer = {}
    for powers in permutations(range(5)):
        alpha = [0] * N
        if raw:
            for i in triple:
                alpha[i] = 2
        for i, power in zip(complement, powers):
            alpha[i] = power
        answer[tuple(alpha)] = sign(powers)
    return answer


def rank_mod(rows, prime=1000003):
    pivots = {}
    for original in rows:
        row = [entry % prime for entry in original]
        for i in range(len(row)):
            entry = row[i]
            if not entry:
                continue
            if i in pivots:
                pivot = pivots[i]
                row = [(a - entry * b) % prime for a, b in zip(row, pivot)]
            else:
                inverse = pow(entry, -1, prime)
                pivots[i] = [a * inverse % prime for a in row]
                break
    return len(pivots)


def coefficient_rows(polynomials):
    support = set().union(*(set(polynomial) for polynomial in polynomials))
    return [[polynomial.get(alpha, 0) for polynomial in polynomials] for alpha in sorted(support)]


def exact_zero_sum_relations(polynomials):
    relations = []
    for i, j in combinations(range(N), 2):
        row = [0] * len(TRIPLES)
        for k in range(N):
            if k in (i, j):
                continue
            triple = tuple(sorted((i, j, k)))
            row[TRIPLES.index(triple)] = (-1) ** (sum(triple) + (k > i) + (k > j))
        result = {}
        for scalar, polynomial in zip(row, polynomials):
            for alpha, coefficient in polynomial.items():
                result[alpha] = result.get(alpha, 0) + scalar * coefficient
        assert all(value == 0 for value in result.values())
        relations.append(row)
    assert rank_mod(relations) == 21
    return relations


def main():
    raw = [coordinate(triple, True) for triple in TRIPLES]
    leading = [coordinate(triple, False) for triple in TRIPLES]
    assert len(CUTS) == 127
    expected_minima = [0, 1, 3, 5, 7, 9, 12]
    for s in range(1, N):
        subset = set(range(s))
        support_minima = [min(sum(alpha[i] for i in subset) for alpha in polynomial)
                          for polynomial in raw]
        formula = [2 * len(subset.intersection(triple))
                   + comb(s - len(subset.intersection(triple)), 2)
                   for triple in TRIPLES]
        assert support_minima == formula
        assert min(formula) == expected_minima[s - 1]
    minima = dict(enumerate(expected_minima, start=1))
    for triple, polynomial in zip(TRIPLES, raw):
        assert len(polynomial) == 120
        assert all(sum(alpha) == 16 and max(alpha) <= 4 for alpha in polynomial)
        assert all(polynomial[tuple(4 - a for a in alpha)] == coefficient
                   for alpha, coefficient in polynomial.items())
        width_sum = 0
        order_sum = 0
        for subset in CUTS:
            support_min = min(sum(alpha[i] for i in subset) for alpha in polynomial)
            centered_max = max(sum(alpha[i] - 2 for i in subset) for alpha in polynomial)
            width_sum += centered_max
            order_sum += support_min - minima[len(subset)]
            assert centered_max + support_min == 2 * len(subset)
        assert (width_sum, order_sum) == (320, 53)
    assert sum(2 * len(s) - minima[len(s)] for s in CUTS) == 373
    relations = exact_zero_sum_relations(raw)
    assert exact_zero_sum_relations(leading) == relations
    assert rank_mod(coefficient_rows(raw)) == 35
    assert rank_mod(coefficient_rows(leading)) == 35
    print("PASS: all 56 raw coordinates: degree 16, separate degree <=4, reciprocal coefficient symmetry.")
    print("PASS: independent support minima for all cut sizes; m=0,1,3,5,7,9,12.")
    print("PASS: every coordinate has width sum 320, normalized-order sum 53; total baseline 373.")
    print("PASS: exact 21-dimensional common relation space and rank-35 witnesses for raw/leading maps.")
    print("CONSEQUENCE: every nonzero degree-one pullback has diagonal multiplicity exactly 10.")


if __name__ == "__main__":
    main()
