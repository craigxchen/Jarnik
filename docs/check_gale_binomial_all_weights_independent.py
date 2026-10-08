#!/usr/bin/env python3
"""Independent factor census, with literal Newton expansions as cross-checks."""

from collections import Counter, defaultdict
from itertools import combinations, combinations_with_replacement

from check_gale_torus_newton_cut_orders_independent import coordinate


LABELS = range(8)
TRIPLES = list(combinations(LABELS, 3))
EDGES = list(combinations(LABELS, 2))
CUTS = [frozenset(s) for k in range(1, 5) for s in combinations(LABELS, k)
        if k != 4 or 0 in s]
BASE = {1: 0, 2: 1, 3: 3, 4: 5}
PRODUCTS = list(combinations_with_replacement(TRIPLES, 2))


def factorization(triples):
    counts = tuple(sum(i in triple for triple in triples) for i in LABELS)
    edge_counts = tuple(sum(i not in triple and j not in triple for triple in triples)
                        for i, j in EDGES)
    return counts, edge_counts


def initial_data(triples, cut):
    counts, edges = factorization(triples)
    powers = [2 * c for c in counts]
    remaining = []
    negative = 0
    for (i, j), multiplicity in zip(EDGES, edges):
        if (i in cut) == (j in cut):
            remaining.append(multiplicity)
        else:
            remaining.append(0)
            powers[j if i in cut else i] += multiplicity
            negative += multiplicity * int(j in cut)
    degree = sum(powers[i] for i in cut)
    degree += sum(m for (i, j), m in zip(EDGES, remaining) if i in cut and j in cut)
    for i in LABELS:
        if i in cut:
            assert powers[i] == 2 * counts[i]
        else:
            maximum = powers[i] + sum(m for edge, m in zip(EDGES, remaining) if i in edge)
            assert maximum == 4 * len(triples) - 2 * counts[i]
    return tuple(powers) + tuple(remaining), (-1) ** negative, degree - len(triples) * BASE[len(cut)]


def polynomial_product(triples):
    result = {(0,) * 8: 1}
    for triple in triples:
        next_result = defaultdict(int)
        for alpha, a in result.items():
            for beta, b in coordinate(triple, True).items():
                next_result[tuple(x + y for x, y in zip(alpha, beta))] += a * b
        result = {alpha: coefficient for alpha, coefficient in next_result.items() if coefficient}
    return result


def literal_orders(left, right, coefficient):
    polynomial = defaultdict(int, polynomial_product(PRODUCTS[left]))
    for alpha, value in polynomial_product(PRODUCTS[right]).items():
        polynomial[alpha] += coefficient * value
    support = [alpha for alpha, value in polynomial.items() if value]
    assert support
    return [min(sum(alpha[i] for i in cut) for alpha in support) - 2 * BASE[len(cut)]
            for cut in CUTS]


def main():
    factorizations = [factorization(product) for product in PRODUCTS]
    assert len(PRODUCTS) == 1596 == len(set(factorizations))
    orders = [[None] * len(CUTS) for _ in PRODUCTS]
    cancellation_cuts = defaultdict(dict)
    for cut_index, cut in enumerate(CUTS):
        groups = defaultdict(list)
        for product_index, product in enumerate(PRODUCTS):
            signature, sign, order = initial_data(product, cut)
            orders[product_index][cut_index] = order
            groups[signature].append((product_index, sign))
        for group in groups.values():
            for (left, a), (right, b) in combinations(group, 2):
                assert orders[left][cut_index] == orders[right][cut_index]
                assert factorizations[left][0] == factorizations[right][0]
                jet = {(i, j): x - y for (i, j), x, y in
                       zip(EDGES, factorizations[left][1], factorizations[right][1])
                       if (i in cut) != (j in cut) and x != y}
                assert jet
                cancellation_cuts[left, right][cut_index] = -a * b
    assert all(sum(row) == 106 for row in orders)
    assert len(cancellation_cuts) == 2100
    assert sum(map(len, cancellation_cuts.values())) == 18480
    histogram = Counter()
    witnesses = {}
    for (left, right), cancellations in cancellation_cuts.items():
        minimum = [min(a, b) for a, b in zip(orders[left], orders[right])]
        for coefficient in (-1, 1):
            predicted = [value + int(cancellations.get(k) == coefficient)
                         for k, value in enumerate(minimum)]
            total = sum(predicted)
            histogram[total] += 1
            witnesses.setdefault(total, (left, right, coefficient, predicted))
    assert histogram == {82: 1260, 86: 1260, 90: 840, 106: 840}
    for total, (left, right, coefficient, predicted) in sorted(witnesses.items()):
        assert literal_orders(left, right, coefficient) == predicted
        print(f"PASS: full Newton expansion at all 127 cuts for order-sum {total}.")
    no_coincidence = next((left, right) for left, right in combinations(range(len(PRODUCTS)), 2)
                          if (left, right) not in cancellation_cuts)
    left, right = no_coincidence
    assert literal_orders(left, right, 2) == [min(a, b) for a, b in zip(orders[left], orders[right])]
    # A degree-four Pasch trade: distinct coordinate monomials become the
    # same polynomial, so the general-degree first-jet lemma needs its
    # explicit distinct-polynomial hypothesis.
    pasch_left = [(0, 1, 2), (0, 3, 4), (1, 3, 5), (2, 4, 5)]
    pasch_right = [(0, 1, 3), (0, 2, 4), (1, 2, 5), (3, 4, 5)]
    assert pasch_left != pasch_right
    assert factorization(pasch_left) == factorization(pasch_right)
    for cut in CUTS:
        assert initial_data(pasch_left, cut) == initial_data(pasch_right, cut)
    print("PASS: independent 1596-product census, 2100 candidate pairs, 18480 coincidences; maximum 106.")
    print("PASS: a noncoinciding coefficient-2 binomial and the degree-four Pasch-trade caveat.")


if __name__ == "__main__":
    main()
