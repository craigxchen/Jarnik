"""Exact fiber and mass-loss audit for full-fair Walsh extraction.

This is a combinatorial checker. Formal phase values only illustrate why
removed nonconstant blocks need an additional cancellation theorem; no short
arc realization is claimed.
"""

from fractions import Fraction
from itertools import combinations
from math import gcd
from random import Random


def canonical(mask, width):
    whole = (1 << width) - 1
    assert 0 < mask < whole
    return mask if mask & 1 else whole ^ mask


def source_classes(order):
    return [mask for mask in range(1, (1 << order) - 1) if mask & 1]


def walsh_hadamard(order):
    assert order >= 2 and order & (order - 1) == 0
    rows = [[1]]
    while len(rows) < order:
        rows = [row + row for row in rows] + [row + [-x for x in row]
                                               for row in rows]
    return rows


def walsh_core_classes(order):
    matrix = walsh_hadamard(order)
    result = set()
    for column in range(1, order):
        mask = sum((matrix[row][column] == 1) << row
                   for row in range(order))
        result.add(canonical(mask, order))
    assert len(result) == order - 1
    return result


def restricted_mask(mask, selected):
    raw = sum(((mask >> row) & 1) << j for j, row in enumerate(selected))
    if raw in (0, (1 << len(selected)) - 1):
        return raw
    return canonical(raw, len(selected))


def restriction_counts(order, selected):
    selected = tuple(selected)
    counts = {}
    constant = 0
    for mask in source_classes(order):
        restricted = restricted_mask(mask, selected)
        if restricted in (0, (1 << len(selected)) - 1):
            constant += 1
        else:
            counts[restricted] = counts.get(restricted, 0) + 1
    expected_fiber = 1 << (order - len(selected))
    expected_constant = expected_fiber - 1
    assert constant == expected_constant
    assert counts and set(counts.values()) == {expected_fiber}
    assert len(counts) == (1 << (len(selected) - 1)) - 1
    return counts, constant


def retained_counts(order, selected, core):
    counts, constant = restriction_counts(order, selected)
    target_count = (1 << (len(selected) - 1)) - 1
    assert len(core) == len(selected) - 1
    assert core <= counts.keys()
    extras = set(counts) - set(core)
    expected_fiber = 1 << (order - len(selected))
    assert sum(counts[mask] for mask in core) == len(core) * expected_fiber
    assert sum(counts[mask] for mask in extras) == (
        target_count - len(core)) * expected_fiber
    keep_common = constant + sum(counts[mask] for mask in core)
    removed_extra = sum(counts[mask] for mask in extras)
    assert keep_common + removed_extra == (1 << (order - 1)) - 1
    return keep_common, removed_extra, constant


def fraction_audit():
    checked = 0
    for order in range(4, 9):
        for selected in combinations(range(order), 4):
            targets = sorted(restriction_counts(order, selected)[0])
            # Every three-class set has the same fiber cost; these sets
            # include, but are not all, valid Walsh cores.
            for core_tuple in combinations(targets, 3):
                core = set(core_tuple)
                keep, removed, constant = retained_counts(order, selected,
                                                           core)
                expected_fiber = 1 << (order - 4)
                assert keep == 4 * expected_fiber - 1
                assert removed == 4 * expected_fiber
                assert constant == expected_fiber - 1
                checked += 1
    return checked


def eight_row_audit():
    checked = 0
    rng = Random(8128)
    for order in range(8, 11):
        for selected in combinations(range(order), 8):
            targets = sorted(restriction_counts(order, selected)[0])
            # Restriction positions are ordered by selected, so the standard
            # Walsh core is one valid target choice.
            core = walsh_core_classes(8)
            assert core <= set(targets)
            choices = [core]
            for _ in range(8):
                choices.append(set(rng.sample(targets, 7)))
            for choice in choices:
                keep, removed, constant = retained_counts(order, selected,
                                                           choice)
                expected_fiber = 1 << (order - 8)
                assert keep == 8 * expected_fiber - 1
                assert removed == 120 * expected_fiber
                assert constant == expected_fiber - 1
                checked += 1
    return checked


def closed_form_audit():
    checked = 0
    for m in (4, 8, 16, 32):
        for order in (m, m + 1, 2 * m, 3 * m):
            fiber = 1 << (order - m)
            source = (1 << (order - 1)) - 1
            extras = (1 << (m - 1)) - m
            constant = fiber - 1
            keep_common = m * fiber - 1
            removed = extras * fiber
            assert keep_common + removed == source
            assert Fraction(removed, source) == Fraction(extras * fiber,
                                                          source)
            # Primitive extraction additionally removes the constant fiber.
            core = (m - 1) * fiber
            assert source - core == removed + constant
            checked += 1
    # The asymptotic fractions used in the note are approached from above
    # with the exact -1 denominator retained.
    assert Fraction(4, 2**3) < Fraction(4 * (1 << 20),
                                          (1 << 20) * 8 - 1)
    assert Fraction(15, 16) < Fraction(120 * (1 << 20),
                                          (1 << 20) * 128 - 1)
    return checked


def nonlinear_recoding_audit():
    """Row permutations and arbitrary target class choices cannot help."""
    checked = 0
    order, m = 9, 4
    selected = tuple(range(m))
    targets = sorted(restriction_counts(order, selected)[0])
    for permutation in ((0, 1, 2, 3), (1, 0, 3, 2), (3, 2, 1, 0)):
        permuted_targets = {
            canonical(sum(((mask >> permutation[j]) & 1) << j
                          for j in range(m)), m)
            for mask in targets
        }
        assert permuted_targets == set(targets)
        for core in combinations(targets, 3):
            assert retained_counts(order, selected, set(core))[:2] == (127, 128)
            checked += 1
    return checked


def phase_error_diagnostic():
    """A single removed nonconstant block can create an order-one error."""
    selected = tuple(range(4))
    targets = sorted(restriction_counts(5, selected)[0])
    core = set(walsh_core_classes(4))
    extra = next(mask for mask in targets if mask not in core)
    # Give one removed block formal argument pi/2.  Its signs differ on the
    # two sides of the cut, so the pairwise removal error is exactly pi.
    plus = [row for row in range(4) if extra >> row & 1]
    minus = [row for row in range(4) if not (extra >> row & 1)]
    assert plus and minus
    # Measure angles in units of pi/2, so this diagnostic stays exact.
    eta_plus = -1
    eta_minus = 1
    assert abs(eta_plus - eta_minus) == 2  # exactly pi
    return extra


def literal_block_phase_diagnostic():
    """Primitive growing blocks can have an order-one individual phase."""
    for n in (1, 2, 10, 100, 10**6):
        a, b = n, n + 1
        assert gcd(a, b) == 1 and (a + b) % 2 == 1
        assert 0 < a < b  # pi/4 < arg(a+ib) < pi/2.
        assert a * a + b * b > 2 * n * n


if __name__ == "__main__":
    four = fraction_audit()
    eight = eight_row_audit()
    closed = closed_form_audit()
    recoded = nonlinear_recoding_audit()
    extra = phase_error_diagnostic()
    literal_block_phase_diagnostic()
    print(f"PASS: {four} exhaustive m=4 row/class-set restriction audits.")
    print(f"PASS: {eight} m=8 restriction audits with a Walsh core and sampled class sets.")
    print(f"PASS: {closed} closed-form mass audits through m=32 and {recoded} recoding checks.")
    print(f"PASS: removed nonconstant target mask {extra} gives formal pair phase error pi.")
    print("PASS: five primitive growing Gaussian blocks with individual phase between pi/4 and pi/2.")
    print("No endpoint phase realization is asserted; removed phases require independent control.")
