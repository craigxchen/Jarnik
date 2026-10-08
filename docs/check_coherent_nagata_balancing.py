#!/usr/bin/env python3
"""Exact coefficient audits for the odd-block balancing identity.

No third-party dependencies. Polynomials in the squarefree z variables
are sparse dictionaries keyed by frozensets of labels; all coefficients
are integers or Fractions. These finite tests audit the uniform proof in
coherent_nagata_balancing.md, not the external generator theorem.
"""
from fractions import Fraction
from itertools import combinations
from math import comb, prod


def clean(poly):
    return {mon: value for mon, value in poly.items() if value}


def add_scaled(target, source, scalar=1):
    for mon, value in source.items():
        target[mon] = target.get(mon, 0) + scalar * value


def multiply(left, right):
    out = {}
    for a, av in left.items():
        for b, bv in right.items():
            assert not a & b, "only disjoint-support products are needed"
            mon = a | b
            out[mon] = out.get(mon, 0) + av * bv
    return clean(out)


def vandermonde(labels, b):
    return prod(b[labels[j]] - b[labels[i]]
                for j in range(len(labels)) for i in range(j))


def odd_determinant(labels, b):
    """Rows 1,b,...,b^s,z,bz,...,b^(s-1)z, columns in labels order."""
    assert len(labels) % 2 == 1 and len(set(labels)) == len(labels)
    s = (len(labels) - 1) // 2
    out = {}
    for positions in combinations(range(len(labels)), s):
        chosen = tuple(labels[p] for p in positions)
        others = tuple(label for p, label in enumerate(labels)
                       if p not in positions)
        sign = (-1) ** (s * (3 * s + 1) // 2 + sum(positions))
        out[frozenset(chosen)] = (sign * vandermonde(chosen, b)
                                 * vandermonde(others, b))
    return clean(out)


def product_of_blocks(blocks, b):
    out = {frozenset(): 1}
    for block in blocks:
        out = multiply(out, odd_determinant(block, b))
    return out


def balancing_terms(block, singleton, b):
    """Return scalar, triple, smaller block for the solved identity (4)."""
    a, *tail = block
    s = len(tail) // 2
    denominator = b[singleton] - b[a]
    assert denominator
    for r, label in enumerate(tail, 1):
        rest = tuple(x for x in tail if x != label)
        numerator = ((-1) ** (r + s - 1)
                     * prod(b[x] - b[a] for x in rest))
        yield (Fraction(numerator, denominator),
               (a, singleton, label), rest)


def check_identity(block, singleton, b):
    target = odd_determinant(block, b)
    actual = {}
    for scalar, triple, rest in balancing_terms(block, singleton, b):
        add_scaled(actual,
                   multiply(odd_determinant(triple, b),
                            odd_determinant(rest, b)), scalar)
    assert clean(actual) == target, (block, singleton)


def recursively_balance(blocks, b):
    """Expand an odd-block product into disjoint-triple/singleton products."""
    oversized = next((i for i, block in enumerate(blocks)
                      if len(block) >= 5), None)
    if oversized is None:
        assert all(len(block) in (1, 3) for block in blocks)
        return [(Fraction(1), blocks)]
    singleton_index = next(i for i, block in enumerate(blocks)
                           if len(block) == 1)
    block = blocks[oversized]
    singleton = blocks[singleton_index][0]
    other = tuple(v for i, v in enumerate(blocks)
                  if i not in (oversized, singleton_index))
    out = []
    for scalar, triple, rest in balancing_terms(block, singleton, b):
        for smaller_scalar, smaller_blocks in recursively_balance(
                other + (triple, rest), b):
            out.append((scalar * smaller_scalar, smaller_blocks))
    return out


def derivative(poly, weights):
    out = {}
    for mon, coefficient in poly.items():
        for i in mon:
            sub = mon - {i}
            out[sub] = out.get(sub, 0) + weights[i] * coefficient
    return clean(out)


def bounded_walk_count(m, d):
    limit = m - 2 * d
    assert limit >= 0
    counts = [1] + [0] * limit
    for _ in range(m):
        following = [0] * (limit + 1)
        for height, value in enumerate(counts):
            if height:
                following[height - 1] += value
            if height < limit:
                following[height + 1] += value
        counts = following
    return counts[limit]


def safe_comb(n, k):
    return comb(n, k) if 0 <= k <= n else 0


def reflected_walk_count(m, d):
    period = m - 2 * d + 2
    return sum(safe_comb(m, d + j * period)
               - safe_comb(m, d - 1 + j * period)
               for j in range(-m - 1, m + 2))


def word_to_chains(m, y_positions):
    """Return the greedy odd-chain decomposition, or None outside the strip."""
    limit = m - 2 * len(y_positions)
    height = 0
    active, inactive = [], []
    for position in range(m):
        if position in y_positions:
            height -= 1
            if height < 0:
                return None
            chain = active.pop()
            inactive.append(chain + [position])
        else:
            height += 1
            if height > limit:
                return None
            if inactive:
                chain = inactive.pop()
                active.append(chain + [position])
            else:
                active.append([position])
    assert height == limit and not inactive and len(active) == limit
    assert sorted(i for chain in active for i in chain) == list(range(m))
    for chain in active:
        assert len(chain) % 2 == 1
        assert all((position in y_positions) == (offset % 2 == 1)
                   for offset, position in enumerate(chain))
    return active


def main():
    identity_count = 0
    for s in range(1, 6):
        count = 2 * s + 2
        fixtures = [
            {i: i * i + 3 * i + 2 for i in range(count)},
            {i: Fraction(2 * i + 1, i + 2) for i in range(count)},
            {i: (-1) ** i * (i + 1) for i in range(count)},
        ]
        for b in fixtures:
            labels = tuple(range(count - 1))
            orders = (labels, labels[::-1], labels[1:] + labels[:1])
            for block in orders:
                check_identity(block, count - 1, b)
                determinant = odd_determinant(block, b)
                assert not derivative(determinant, {i: 1 for i in b})
                assert not derivative(determinant, b)
                identity_count += 1

    b = {i: Fraction(i * i + 2 * i + 3, i + 2) for i in range(12)}
    sizes_to_check = ((9, 1, 1, 1), (7, 3, 1, 1),
                      (5, 5, 1, 1), (5, 3, 3, 1), (3, 3, 3, 3))
    term_counts = []
    for sizes in sizes_to_check:
        start = 0
        blocks = []
        for size in sizes:
            blocks.append(tuple(range(start, start + size)))
            start += size
        assert start == 12
        expanded = recursively_balance(tuple(blocks), b)
        actual = {}
        for scalar, triple_blocks in expanded:
            assert sorted(map(len, triple_blocks)) == [3] * 4
            add_scaled(actual, product_of_blocks(triple_blocks, b), scalar)
        target = product_of_blocks(blocks, b)
        assert clean(actual) == target
        term_counts.append((sizes, len(expanded)))

    # The below-one-third version leaves unused singleton labels.
    blocks = ((0, 1, 2, 3, 4), (5,), (6,), (7,))
    expanded = recursively_balance(blocks, b)
    actual = {}
    for scalar, balanced in expanded:
        assert sorted(map(len, balanced)) == [1, 1, 3, 3]
        add_scaled(actual, product_of_blocks(balanced, b), scalar)
    assert clean(actual) == product_of_blocks(blocks, b)

    for m in range(1, 13):
        for d in range(m // 2 + 1):
            count = sum(word_to_chains(m, frozenset(y)) is not None
                        for y in combinations(range(m), d))
            assert count == bounded_walk_count(m, d)
            assert count == reflected_walk_count(m, d)
    ranks = []
    for k in range(1, 9):
        m, d = 6 * k, 2 * k
        expected = comb(m, d) - 2 * comb(m, d - 1) + comb(m, d - 2)
        assert bounded_walk_count(m, d) == reflected_walk_count(m, d) == expected
        ranks.append(expected)

    print(f"PASS: {identity_count} exact determinant identities and invariance checks")
    print("PASS: twelve-label recursive balancing", term_counts)
    print("PASS: below-one-third balancing with unused labels")
    print("PASS: odd-chain words and reflection formula for m<=12")
    print("PASS: rectangular terminal ranks for k=1..8", ranks)


if __name__ == "__main__":
    main()
