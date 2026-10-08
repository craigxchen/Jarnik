"""Elementary affine cubes and exact Gaussian whole-cut restriction.

The literal Gaussian fixtures do not have controlled arguments. The
effective endpoint exclusion is the previously proved arithmetic input.
"""

from fractions import Fraction
from functools import reduce
from itertools import combinations
from random import Random

from check_five_row_gradient_cofactor_stress import (
    conj, gaussian_gcd, mul, norm, power,
)


def density_condition(size, order, dimension):
    exponent = 2 ** (dimension - 1)
    return size ** exponent >= exponent * 2 ** exponent * order ** (exponent - 1)


def dense_affine_flat(points, order, dimension):
    assert order & (order - 1) == 0 and dimension >= 1
    assert points and all(0 <= x < order for x in points)
    assert density_condition(len(points), order, dimension)
    current, linear = set(points), {0}
    lower_density = Fraction(len(points), order)
    for i in range(dimension):
        assert len(linear) == 2 ** i
        assert all(x ^ u in current for x in current for u in linear)
        assert Fraction(len(current), order) >= lower_density
        assert len(current) >= 2 ** (i + 1)
        choices = [(len(current & {x ^ h for x in current}), h)
                   for h in range(order) if h not in linear]
        assert sum(count for count, _ in choices) == len(current) ** 2 - len(linear) * len(current)
        count, direction = max(choices)
        assert 2 * order * count >= len(current) ** 2
        current &= {x ^ direction for x in current}
        linear |= {u ^ direction for u in linear}
        lower_density = lower_density ** 2 / 2
    assert current and len(linear) == 2 ** dimension
    base = min(current)
    flat = {base ^ u for u in linear}
    assert len(flat) == 2 ** dimension and flat <= points
    return flat


def affine_cube_fixtures():
    checked = 0
    # Exhaust the small density range covered by the exact dimension-two
    # bound, rather than searching arbitrary subsets indefinitely.
    for size in range(12, 17):
        for subset in combinations(range(16), size):
            assert density_condition(size, 16, 2)
            dense_affine_flat(set(subset), 16, 2)
            checked += 1
    rng = Random(278)
    for order, dimension, size in ((32, 2, 20), (64, 2, 32),
                                   (128, 3, 112), (256, 3, 192)):
        for _ in range(4):
            points = set(rng.sample(range(order), size))
            dense_affine_flat(points, order, dimension)
            checked += 1
    return checked


def character(v, x):
    return -1 if bin(v & x).count("1") % 2 else 1


def canonical_mask(mask, count):
    whole = (1 << count) - 1
    assert mask not in (0, whole)
    return mask if mask & 1 else whole ^ mask


def walsh_cut_classes(labels, ambient_order):
    count = len(labels)
    whole = (1 << count) - 1
    classes = set()
    for v in range(1, ambient_order):
        mask = sum((character(v, x) == 1) << j for j, x in enumerate(labels))
        if mask not in (0, whole):
            classes.add(canonical_mask(mask, count))
    return classes


def affine_planes(order):
    flats = set()
    for first, second in combinations(range(1, order), 2):
        linear = {0, first, second, first ^ second}
        for base in range(order):
            flats.add(tuple(sorted(base ^ x for x in linear)))
    return sorted(flats)


def threshold_factorization(common, primes, allocations, widths, labels):
    content, blocks = common, {}
    whole = (1 << len(labels)) - 1
    for prime, column, width in zip(primes, allocations, widths):
        for level in range(1, width + 1):
            raw_mask = sum((column[x] >= level) << j for j, x in enumerate(labels))
            if raw_mask == 0:
                content = mul(content, conj(prime))
            elif raw_mask == whole:
                content = mul(content, prime)
            else:
                mask = canonical_mask(raw_mask, len(labels))
                oriented_prime = prime if raw_mask == mask else conj(prime)
                blocks[mask] = mul(blocks.get(mask, (1, 0)), oriented_prime)
    return content, blocks


def merge_restricted_blocks(common, blocks, old_labels, labels):
    content, restricted = common, {}
    old_indices = {x: j for j, x in enumerate(old_labels)}
    whole = (1 << len(labels)) - 1
    for old_mask, block in blocks.items():
        raw_mask = sum(bool(old_mask & (1 << old_indices[x])) << j
                       for j, x in enumerate(labels))
        if raw_mask == 0:
            content = mul(content, conj(block))
        elif raw_mask == whole:
            content = mul(content, block)
        else:
            mask = canonical_mask(raw_mask, len(labels))
            oriented = block if raw_mask == mask else conj(block)
            restricted[mask] = mul(restricted.get(mask, (1, 0)), oriented)
    return content, restricted


def gaussian_restriction_fixture():
    order = 8
    primes = [(2, 1), (3, 2), (4, 1), (5, 2), (6, 1), (5, 4), (7, 2)]
    common = (7, 11)
    allocations, widths = [], []
    for v in range(1, order):
        middle = {x for x in range(order) if character(v, x) == 1}
        if v == 1:
            chain = [{0}, middle, middle | {1}]
        elif v == 2:
            chain = [{0, 1}, middle, middle | {2, 3}]
        else:
            chain = [middle] * (1 + (v % 2))
        assert all(first <= second for first, second in zip(chain, chain[1:]))
        column = [sum(x in cut for cut in chain) for x in range(order)]
        assert min(column) == 0 and max(column) == len(chain)
        allocations.append(column)
        widths.append(len(chain))
    units = [power((0, 1), (3 * x + 1) % 4) for x in range(order)]
    rows = []
    for x in range(order):
        row = mul(common, units[x])
        for prime, column, width in zip(primes, allocations, widths):
            row = mul(row, power(prime, column[x]))
            row = mul(row, power(conj(prime), width - column[x]))
        rows.append(row)
    assert len(set(rows)) == order
    source_norm = norm(rows[0])
    assert all(norm(row) == source_norm for row in rows)
    old_labels = tuple(range(order))
    old_content, old_blocks = threshold_factorization(common, primes, allocations, widths, old_labels)
    assert old_content == common
    core_classes = walsh_cut_classes(old_labels, order)
    assert len(core_classes) == order - 1 and core_classes <= old_blocks.keys()
    extra_count = len(old_blocks) - len(core_classes)
    assert extra_count == 4
    checked, content_growth, orientation_changes = 0, 0, 0
    for labels in affine_planes(order):
        content, blocks = merge_restricted_blocks(old_content, old_blocks, old_labels, labels)
        direct_content, direct_blocks = threshold_factorization(common, primes, allocations, widths, labels)
        assert content == direct_content and blocks == direct_blocks
        assert all(norm(block) >= 5 and norm(gaussian_gcd(block, conj(block))) == 1
                   for block in blocks.values())
        restricted_core = walsh_cut_classes(labels, order)
        assert len(restricted_core) == len(labels) - 1 == 3
        assert restricted_core <= blocks.keys()
        assert all(sum((1 if mask & (1 << i) else -1)
                       * (1 if mask & (1 << j) else -1)
                       for mask in restricted_core) == -1
                   for i, j in combinations(range(len(labels)), 2))
        assert len(blocks) - len(restricted_core) <= extra_count
        for j, x in enumerate(labels):
            reconstructed = mul(content, units[x])
            for mask, block in blocks.items():
                reconstructed = mul(reconstructed, block if mask & (1 << j) else conj(block))
            assert reconstructed == rows[x]
            assert norm(reconstructed) == source_norm
        actual_gcd = reduce(gaussian_gcd, (rows[x] for x in labels))
        assert norm(actual_gcd) == norm(content)
        assert norm(content) * reduce(lambda a, b: a * b, (norm(g) for g in blocks.values()), 1) == source_norm
        content_growth += norm(content) > norm(common)
        for mask in old_blocks:
            raw_mask = sum(bool(mask & (1 << x)) << j for j, x in enumerate(labels))
            orientation_changes += raw_mask not in (0, 15) and not (raw_mask & 1)
        # Same-prime raw threshold chains never contribute complementary
        # nonconstant cuts on the retained rows.
        for column, width in zip(allocations, widths):
            masks = [sum((column[x] >= level) << j for j, x in enumerate(labels))
                     for level in range(1, width + 1)]
            assert not any(first not in (0, 15) and second == (15 ^ first)
                           for first in masks for second in masks)
        checked += 1
    assert checked == 14 and content_growth and orientation_changes
    # Arbitrary formal block factors do not have this property.
    pi = (2, 1)
    formal_columns = {0b011: pi, 0b001: conj(pi)}
    formal_content, formal_blocks = merge_restricted_blocks((1, 0), formal_columns, (0, 1, 2), (0, 2))
    assert formal_content == (1, 0) and list(formal_blocks.values()) == [(5, 0)]
    assert norm(gaussian_gcd((5, 0), (5, 0))) == 25
    return checked, content_growth, orientation_changes


def full_fair_sparse_simplex():
    checked = 0
    previous_density = Fraction(1)
    for m in range(3, 9):
        order = 2 ** (m - 1)
        labels = (0,) + tuple(1 << j for j in range(m - 1))
        classes = walsh_cut_classes(labels, order)
        assert len(classes) == 2 ** (m - 1) - 1
        assert classes == {mask for mask in range(1, (1 << m) - 1) if mask & 1}
        assert all(a ^ b ^ c ^ d != 0 for a, b, c, d in combinations(labels, 4))
        density = Fraction(m, order)
        assert density < previous_density
        previous_density = density
        checked += 1
    return checked


if __name__ == "__main__":
    cubes = affine_cube_fixtures()
    planes, growth, flips = gaussian_restriction_fixture()
    sparse = full_fair_sparse_simplex()
    print(f"PASS: {cubes} exact dense affine-flat constructions and averaging identities.")
    print(f"PASS: {planes} Gaussian affine-plane restrictions; {growth} enlarged common contents and {flips} orientation reversals.")
    print("PASS: nested-prime noncancellation, literal row units, unchanged radius, and complete nonunit Walsh cores.")
    print(f"PASS: {sparse} sparse full-fair simplices with no affine two-flat.")
    print("Gaussian arguments are uncontrolled; effective endpoint exclusion is a cited proof input.")
