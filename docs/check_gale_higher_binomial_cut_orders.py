#!/usr/bin/env python3
"""Exact all-cut scan of cubic binomials in eight-point Gale coordinates."""

from collections import Counter, defaultdict
from itertools import combinations, combinations_with_replacement


LABELS = tuple(range(8))
TRIPLES = list(combinations(LABELS, 3))
CUTS = [
    cut
    for size in range(1, 5)
    for cut in combinations(LABELS, size)
    if size < 4 or 0 in cut
]
BASELINES = (0, 1, 3, 5)


def incidence(product):
    answer = [0] * 8
    for triple in product:
        for label in triple:
            answer[label] += 1
    return tuple(answer)


def edge_multiplicities(product):
    return tuple(
        sum(i not in triple and j not in triple for triple in product)
        for i, j in combinations(LABELS, 2)
    )


def factor_signature(product):
    """Unique-factorization data for the full phase polynomial product."""
    return incidence(product) + edge_multiplicities(product)


def coordinate_raw_order(triple, cut):
    intersection = len(set(triple) & set(cut))
    small_complement = len(cut) - intersection
    return 2 * intersection + small_complement * (small_complement - 1) // 2


def normalized_orders(product):
    return tuple(
        sum(coordinate_raw_order(triple, cut) for triple in product)
        - 3 * BASELINES[len(cut) - 1]
        for cut in CUTS
    )


def coordinate_initial_signature(triple, cut):
    """Sign, variable powers, and same-side difference powers at a cut."""
    cut = set(cut)
    complement = set(LABELS) - set(triple)
    variable_powers = [0] * 8
    difference_powers = [0] * 28
    sign = 1
    for label in triple:
        variable_powers[label] += 2
    for position, (i, j) in enumerate(combinations(LABELS, 2)):
        if i not in complement or j not in complement:
            continue
        if (i in cut) == (j in cut):
            difference_powers[position] += 1
        else:
            outside = j if i in cut else i
            variable_powers[outside] += 1
            # The factor is z_j-z_i.  Its leading unit is negative exactly
            # when the lower-index endpoint is outside the cut.
            if i not in cut:
                sign = -sign
    return sign, tuple(variable_powers), tuple(difference_powers)


def product_initial_signature(product, cut):
    sign = 1
    variable_powers = [0] * 8
    difference_powers = [0] * 28
    for triple in product:
        triple_sign, variables, differences = coordinate_initial_signature(triple, cut)
        sign *= triple_sign
        variable_powers = [x + y for x, y in zip(variable_powers, variables)]
        difference_powers = [x + y for x, y in zip(difference_powers, differences)]
    return sign, tuple(variable_powers), tuple(difference_powers)


def all_products():
    return [
        tuple(TRIPLES[index] for index in indices)
        for indices in combinations_with_replacement(range(len(TRIPLES)), 3)
    ]


def scan_canonical_weight(product_group):
    orders = [normalized_orders(product) for product in product_group]
    initials = [
        [product_initial_signature(product, cut) for cut in CUTS]
        for product in product_group
    ]
    assert all(sum(order) == 159 for order in orders)

    maximum = -1
    histogram = Counter()
    equality_count = 0
    for left, right in combinations(range(len(product_group)), 2):
        left_orders = orders[left]
        right_orders = orders[right]
        absolute_difference = sum(
            abs(a - b) for a, b in zip(left_orders, right_orders)
        )
        assert absolute_difference % 2 == 0
        common_minimum = 159 - absolute_difference // 2

        cancellation_for_plus = 0
        cancellation_for_minus = 0
        for position, (a, b) in enumerate(zip(left_orders, right_orders)):
            if a != b:
                continue
            left_initial = initials[left][position]
            right_initial = initials[right][position]
            if left_initial[1:] != right_initial[1:]:
                continue
            # Distinct full products with equal cut initials differ in a
            # cross-edge multiplicity.  Their first logarithmic jets are
            # therefore distinct, so cancellation raises the order by one.
            assert factor_signature(product_group[left]) != factor_signature(
                product_group[right]
            )
            if left_initial[0] == right_initial[0]:
                cancellation_for_minus += 1
            else:
                cancellation_for_plus += 1

        plus_order = common_minimum + cancellation_for_plus
        minus_order = common_minimum + cancellation_for_minus
        histogram[plus_order] += 1
        histogram[minus_order] += 1
        maximum = max(maximum, plus_order, minus_order)
        equality_count += (plus_order == 159) + (minus_order == 159)
    return maximum, histogram, equality_count


def main():
    assert len(CUTS) == 127
    products = all_products()
    assert len(products) == 30_856

    # No two cubic coordinate products are the same phase polynomial up to
    # their fixed orientation sign.  The first factor trades occur later.
    signatures = {factor_signature(product) for product in products}
    assert len(signatures) == len(products)

    weight_groups = defaultdict(list)
    for product in products:
        weight_groups[incidence(product)].append(product)
    assert len(weight_groups) == 5_328

    weight_types = defaultdict(list)
    for weight, group in weight_groups.items():
        weight_types[tuple(sorted(weight, reverse=True))].append((weight, group))
    expected_sizes = {
        (3, 3, 3, 0, 0, 0, 0, 0): 1,
        (3, 3, 2, 1, 0, 0, 0, 0): 1,
        (3, 3, 1, 1, 1, 0, 0, 0): 1,
        (3, 2, 2, 2, 0, 0, 0, 0): 1,
        (3, 2, 2, 1, 1, 0, 0, 0): 3,
        (3, 2, 1, 1, 1, 1, 0, 0): 6,
        (3, 1, 1, 1, 1, 1, 1, 0): 15,
        (2, 2, 2, 2, 1, 0, 0, 0): 6,
        (2, 2, 2, 1, 1, 1, 0, 0): 16,
        (2, 2, 1, 1, 1, 1, 1, 0): 40,
        (2, 1, 1, 1, 1, 1, 1, 1): 105,
    }
    assert set(weight_types) == set(expected_sizes)

    global_maximum = -1
    signed_pairs = 0
    equality_count = 0
    total_histogram = Counter()
    for weight_type, expected_size in expected_sizes.items():
        representatives = weight_types[weight_type]
        canonical_group = next(
            group for weight, group in representatives if weight == weight_type
        )
        assert len(canonical_group) == expected_size
        assert all(len(group) == expected_size for _, group in representatives)
        if expected_size == 1:
            continue
        maximum, histogram, equalities = scan_canonical_weight(canonical_group)
        global_maximum = max(global_maximum, maximum)
        total_histogram.update(histogram)
        signed_pairs += expected_size * (expected_size - 1)
        equality_count += equalities

    assert signed_pairs == 12_996
    assert sum(total_histogram.values()) == signed_pairs
    assert global_maximum == 159
    assert equality_count == 264
    assert max(total_histogram) == 159
    print(
        "PASS: all canonical cubic Gale binomials have total clean-cut "
        "order at most 159."
    )


if __name__ == "__main__":
    main()
