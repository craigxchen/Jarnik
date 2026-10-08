"""Exact factor-based cut orders for every eight-point quadratic binomial.

No polynomial sampling is used. Unique factorization gives leading-form
equality; independence of the cross ratios gives exactly one extra order
when two distinct products have cancelling leading forms.
"""

from collections import defaultdict
from itertools import combinations, combinations_with_replacement


N = 8
PAIRS = list(combinations(range(N), 2))
TRIPLES = list(combinations(range(N), 3))
CUTS = [set(cut) for size in (1, 2, 3) for cut in combinations(range(N), size)]
CUTS += [set(cut) for cut in combinations(range(N), 4) if 0 in cut]
BASELINES = {1: 0, 2: 1, 3: 3, 4: 5}


def raw_factors(triple):
    inside = set(triple)
    return (tuple(2 * int(i in inside) for i in range(N)),
            tuple(int(i not in inside and j not in inside) for i, j in PAIRS))


RAW = [raw_factors(triple) for triple in TRIPLES]
PRODUCTS = list(combinations_with_replacement(range(len(TRIPLES)), 2))
FACTORS = [(tuple(a + b for a, b in zip(RAW[i][0], RAW[j][0])),
            tuple(a + b for a, b in zip(RAW[i][1], RAW[j][1])))
           for i, j in PRODUCTS]


def cut_data(factors, cut):
    powers, edges = factors
    leading_powers = list(powers)
    retained_edges = list(edges)
    order = sum(powers[i] for i in cut)
    sign = 1
    for k, (i, j) in enumerate(PAIRS):
        multiplicity = edges[k]
        if i in cut and j in cut:
            order += multiplicity
        elif (i in cut) != (j in cut):
            retained_edges[k] = 0
            outside = j if i in cut else i
            leading_powers[outside] += multiplicity
            if j in cut and multiplicity % 2:
                sign = -sign
    primitive_order = order - 2 * BASELINES[len(cut)]
    assert primitive_order >= 0
    signature = tuple(leading_powers) + tuple(retained_edges)
    return primitive_order, sign, signature


def main():
    assert len(PRODUCTS) == 1596 and len(CUTS) == 127
    assert len(set(FACTORS)) == len(FACTORS), 'distinct products must be distinct polynomials'
    valuations = [[] for _ in PRODUCTS]
    cancellations = defaultdict(lambda: [0, 0])  # coefficients +1, -1
    total_coincident_pairs = 0
    for cut in CUTS:
        groups = defaultdict(list)
        for index, factors in enumerate(FACTORS):
            order, sign, signature = cut_data(factors, cut)
            valuations[index].append(order)
            groups[(order, signature)].append((index, sign))
        for group in groups.values():
            for (left, left_sign), (right, right_sign) in combinations(group, 2):
                # Distinct cross-factor multiplicities give a nonzero first jet.
                assert any(a != b and ((i in cut) != (j in cut))
                           for (i, j), a, b in zip(PAIRS, FACTORS[left][1], FACTORS[right][1]))
                coefficient = -left_sign * right_sign
                cancellations[left, right][0 if coefficient == 1 else 1] += 1
                total_coincident_pairs += 1
    assert all(sum(row) == 106 for row in valuations)

    best = -1
    witnesses = []
    histogram = defaultdict(int)
    type_counts = defaultdict(int)
    for (left, right), counts in cancellations.items():
        assert FACTORS[left][0] == FACTORS[right][0]
        loss_twice = sum(abs(a - b) for a, b in zip(valuations[left], valuations[right]))
        assert loss_twice % 2 == 0
        loss = loss_twice // 2
        repeated_labels = sum(power == 4 for power in FACTORS[left][0])
        if repeated_labels == 1:
            assert loss == 16 and sorted(counts) == [0, 16]
        else:
            assert repeated_labels == 0
            assert loss == 24 and sorted(counts) == [0, 4]
        type_counts[repeated_labels] += 1
        for coefficient, count in zip((1, -1), counts):
            total_order = 106 - loss + count
            histogram[total_order] += 1
            if total_order > best:
                best = total_order
                witnesses = [(left, right, coefficient, loss, count)]
            elif total_order == best and len(witnesses) < 5:
                witnesses.append((left, right, coefficient, loss, count))
    # Every pair omitted from cancellations has no increase over its minimum
    # order, whose total is <=106. Coefficients other than +/-1 cannot cancel.
    print('Products:', len(PRODUCTS), '; cuts:', len(CUTS))
    print('Pairs with a possible leading cancellation:', len(cancellations))
    print('Total leading-form coincidences:', total_coincident_pairs)
    print('Maximum order sum among these signed binomials:', best)
    for left, right, coefficient, loss, count in witnesses:
        print('Witness:', tuple(TRIPLES[i] for i in PRODUCTS[left]), coefficient,
              tuple(TRIPLES[i] for i in PRODUCTS[right]), 'loss', loss, 'gain', count)
    print('Order histogram:', dict(sorted(histogram.items())))
    assert best == 106
    assert dict(type_counts) == {1: 840, 0: 1260}


if __name__ == '__main__':
    main()
