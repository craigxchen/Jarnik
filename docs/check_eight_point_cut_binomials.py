#!/usr/bin/env python3
"""Bounded exact Newton scan for degree-two eight-point Plucker binomials.

The pullback is computed in the phase torus, with

    R_I = (prod_(i in I) z_i^2) prod_(j<k, j,k notin I) (z_k-z_j).

All arithmetic is integral.  For a cut S of minority size s, the primitive
order of a nonzero pullback F is

    min_(alpha in supp(F)) sum_(i in S) alpha_i - 2*m_s,

where m_s=(0,1,3,5).  The script scans the one representative weight block
of sizes 3 and 10, for every signed 2- and 3-term combination.  The size-1
blocks are monomials and are checked separately.  Relabeling carries these
representatives to all blocks of the same type.
"""

from collections import defaultdict
from itertools import combinations, permutations, product


TRIPLES = list(combinations(range(8), 3))
PERMUTATIONS = list(permutations(range(5)))
M_S = {1: 0, 2: 1, 3: 3, 4: 5}


def permutation_sign(values):
    return -1 if sum(values[i] > values[j]
                    for i in range(len(values))
                    for j in range(i + 1, len(values))) % 2 else 1


def vandermonde_support(nodes):
    """Return the signed exponent support of a five-node Vandermonde."""
    result = defaultdict(int)
    for permutation in PERMUTATIONS:
        exponent = [0] * 8
        for node, power in zip(nodes, permutation):
            exponent[node] = power
        result[tuple(exponent)] += permutation_sign(permutation)
    return {exponent: coefficient for exponent, coefficient in result.items()
            if coefficient}


def multiply(left, right):
    result = defaultdict(int)
    for exponent_left, coefficient_left in left.items():
        for exponent_right, coefficient_right in right.items():
            exponent = tuple(a + b for a, b in zip(exponent_left, exponent_right))
            result[exponent] += coefficient_left * coefficient_right
    return {exponent: coefficient for exponent, coefficient in result.items()
            if coefficient}


def add_scaled(target, source, scale):
    for exponent, coefficient in source.items():
        target[exponent] += scale * coefficient
        if not target[exponent]:
            del target[exponent]


def coordinate(triple):
    complement = tuple(i for i in range(8) if i not in triple)
    support = vandermonde_support(complement)
    shift = [0] * 8
    for i in triple:
        shift[i] = 2
    return {tuple(a + b for a, b in zip(exponent, shift)): coefficient
            for exponent, coefficient in support.items()}


def monomial(left, right, coordinates):
    return multiply(coordinates[left], coordinates[right])


def weight(left, right):
    return tuple((i in left) + (i in right) for i in range(8))


def blocks():
    groups = defaultdict(list)
    for position, left in enumerate(TRIPLES):
        for right in TRIPLES[position:]:
            groups[weight(left, right)].append((left, right))
    return groups


def cuts():
    """All 127 unoriented cuts, represented by masks containing label 0."""
    result = []
    for mask in range(1, 256, 2):
        if mask == 255:
            continue
        side = {i for i in range(8) if mask >> i & 1}
        if len(side) > 4:
            side = set(range(8)) - side
        result.append(tuple(sorted(side)))
    assert len(result) == 127
    return result


def order(support, side):
    s = len(side)
    assert 1 <= s <= 4
    return min(sum(exponent[i] for i in side) for exponent in support) - 2 * M_S[s]


def coordinate_order(triple, side):
    """The direct Vandermonde minimum, after local normalization."""
    s = len(side)
    a = len(set(triple) & set(side))
    return 2 * a + (s - a) * (s - a - 1) // 2 - M_S[s]


def representative_groups(groups):
    by_size = defaultdict(list)
    for block in groups.values():
        by_size[len(block)].append(block)
    # r=1 and r=0 representatives from the conic identities in the barrier note.
    r1 = next(block for block in by_size[3]
              if set(block) == {
                  ((0, 1, 2), (0, 3, 4)),
                  ((0, 1, 3), (0, 2, 4)),
                  ((0, 1, 4), (0, 2, 3)),
              })
    r0 = next(block for block in by_size[10]
              if weight(*block[0]) == (1, 1, 1, 1, 1, 1, 0, 0))
    return r1, r0


def scan_block(block, coordinates, all_cuts):
    terms = [monomial(*pair, coordinates) for pair in block]
    # Cache only initial Newton forms.  Full support is needed only when the
    # initial forms cancel; this keeps the 10-term block bounded.
    initials = []
    for support in terms:
        by_cut = []
        for side in all_cuts:
            minimum = min(sum(exponent[i] for i in side)
                          for exponent in support)
            initial = {}
            for exponent, coefficient in support.items():
                if sum(exponent[i] for i in side) == minimum:
                    initial[exponent] = coefficient
            by_cut.append((minimum, initial))
        initials.append(by_cut)

    def combination_order(selected, coefficients, cut_index, fallback):
        data = [(initials[position][cut_index][0],
                 initials[position][cut_index][1], coefficient)
                for position, coefficient in zip(selected, coefficients)]
        minimum = min(item[0] for item in data)
        leading = defaultdict(int)
        for level, initial, coefficient in data:
            if level == minimum:
                add_scaled(leading, initial, coefficient)
        if leading:
            return minimum - 2 * M_S[len(all_cuts[cut_index])]
        return order(fallback, all_cuts[cut_index])

    results = []
    for size in (2, 3):
        for selected in combinations(range(len(block)), size):
            for signs in product((1, -1), repeat=size - 1):
                coefficients = (1,) + signs
                support = defaultdict(int)
                for position, coefficient in zip(selected, coefficients):
                    add_scaled(support, terms[position], coefficient)
                if not support:
                    continue
                support = dict(support)
                total = sum(combination_order(selected, coefficients, cut_index,
                                              support)
                            for cut_index in range(len(all_cuts)))
                results.append((total, selected, coefficients, len(support)))
    return results


def main():
    groups = blocks()
    histogram = defaultdict(int)
    for block in groups.values():
        histogram[len(block)] += 1
    assert dict(histogram) == {1: 476, 3: 280, 10: 28}
    all_cuts = cuts()
    coordinates = {triple: coordinate(triple) for triple in TRIPLES}

    # Every degree-two coordinate product has order budget 106.  Check the
    # support formula on three products and the direct minima on all 1596.
    for pair in [(TRIPLES[0], TRIPLES[0]),
                 (TRIPLES[0], TRIPLES[1]),
                 (TRIPLES[0], TRIPLES[-1])]:
        support = monomial(*pair, coordinates)
        assert sum(order(support, side) for side in all_cuts) == 106
    for left in TRIPLES:
        for right in TRIPLES:
            assert sum(coordinate_order(left, side) +
                       coordinate_order(right, side) for side in all_cuts) == 106

    r1, r0 = representative_groups(groups)
    scans = {3: scan_block(r1, coordinates, all_cuts),
             10: scan_block(r0, coordinates, all_cuts)}
    for size, results in scans.items():
        assert results
        maximum = max(result[0] for result in results)
        counts = defaultdict(int)
        for total, _, _, support_size in results:
            counts[total] += 1
        print(f"size {size}: combinations={len(results)}, max_total={maximum}, "
              f"totals={dict(sorted(counts.items()))}")
        assert maximum == 106
    print("PASS: exact Newton scan finds no nonzero signed 2/3-term gain above 106.")


if __name__ == "__main__":
    main()
