"""Reproducible finite audit of nonlinear Walsh coset gains (M=32,64).

The gain of a coset/restriction is ``2*k-h``, where ``h=|V|`` and
``k`` is the number of rows whose assigned label has that restriction.
All searches below are exhaustive over the listed subspaces.  The M=64
assignment experiments are only fixed-seed examples, not a theorem about
arbitrary assignments.
"""

from collections import Counter
from fractions import Fraction
from itertools import combinations
import random


def parity(n):
    return bin(n).count("1") & 1


def subspaces(order, dimension):
    """Enumerate subspaces of F_2^log2(order), as sorted element tuples."""
    result = set()
    for basis in combinations(range(1, order), dimension):
        elements = {0}
        for vector in basis:
            elements |= {x ^ vector for x in tuple(elements)}
        if len(elements) == 1 << dimension:
            result.add(tuple(sorted(elements)))
    return sorted(result)


def coset_stats(labels, order, directions):
    dimension = len(directions[0]).bit_length() - 1
    h = 1 << dimension
    restrictions_count = h - 1
    best_gains = Counter()
    positive_vectors = Counter()
    positive_cosets = 0
    positive_triples = 0
    total_cosets = 0

    for vector_space in directions:
        patterns = {
            label: tuple(parity(label & v) for v in vector_space)
            for label in range(1, order)
        }
        zero = (0,) * h
        restrictions = sorted(set(patterns.values()) - {zero})
        representatives = {
            min(x ^ v for v in vector_space) for x in range(order)
        }
        for representative in representatives:
            counts = [
                sum(
                    patterns[labels[representative ^ v]] == restriction
                    for v in vector_space
                )
                for restriction in restrictions
            ]
            maximum = max(counts)
            gain = 2 * maximum - h
            best_gains[gain] += 1
            total_cosets += 1
            positive_triples += sum(2 * count > h for count in counts)
            if gain > 0:
                positive_cosets += 1
                positive_vectors[tuple(sorted(counts, reverse=True))] += 1

    return {
        "directions": len(directions),
        "cosets": total_cosets,
        "positive_cosets": positive_cosets,
        "positive_triples": positive_triples,
        "max_gain": max(best_gains),
        "mean_best": Fraction(
            sum(gain * count for gain, count in best_gains.items()),
            total_cosets,
        ),
        "best_gains": dict(sorted(best_gains.items())),
        "positive_vectors": dict(sorted(positive_vectors.items())),
        "nonuniform_positive": all(
            len(set(vector)) != 1 for vector in positive_vectors
        ),
        "restriction_count": restrictions_count,
    }


def random_assignment(seed):
    rng = random.Random(seed)
    labels = list(range(1, 64)) + [1]
    rng.shuffle(labels)
    return labels


def affine_assignment(seed):
    """Full-rank affine map, with its zero output replaced by duplicate 1."""
    rng = random.Random(seed)
    columns = []
    span = {0}
    while len(columns) < 6:
        candidate = rng.randrange(1, 64)
        if candidate not in span:
            columns.append(candidate)
            span |= {x ^ candidate for x in tuple(span)}
    offset = rng.randrange(64)
    labels = []
    for row in range(64):
        value = offset
        for bit, column in enumerate(columns):
            if (row >> bit) & 1:
                value ^= column
        labels.append(1 if value == 0 else value)
    assert len(span) == 64
    return labels, columns, offset


def controlled_swaps(labels, seed, count):
    rng = random.Random(seed)
    labels = list(labels)
    rows = rng.sample(range(64), 2 * count)
    for i in range(count):
        left, right = rows[2 * i:2 * i + 2]
        labels[left], labels[right] = labels[right], labels[left]
    return labels


def winning_restriction_balance(labels, directions):
    """Count directions where positive-coset winners are exactly balanced.

    In these dimensions a positive coset has a unique winning restriction,
    since two winning restrictions would require more than h rows.
    This is a direct tally, not an LP feasibility assertion.
    """
    dimension = len(directions[0]).bit_length() - 1
    h = 1 << dimension
    q = h - 1
    balanced = 0
    for vector_space in directions:
        patterns = {
            label: tuple(parity(label & v) for v in vector_space)
            for label in range(1, 64)
        }
        zero = (0,) * h
        restrictions = sorted(set(patterns.values()) - {zero})
        representatives = {
            min(x ^ v for v in vector_space) for x in range(64)
        }
        winners = [0] * q
        for representative in representatives:
            counts = [
                sum(
                    patterns[labels[representative ^ v]] == restriction
                    for v in vector_space
                )
                for restriction in restrictions
            ]
            if 2 * max(counts) > h:
                assert counts.count(max(counts)) == 1
                winners[counts.index(max(counts))] += 1
        if sum(winners) and len(set(winners)) == 1:
            balanced += 1
    return balanced


def main():
    m32 = (
        28, 21, 11, 30, 24, 23, 13, 3,
        12, 7, 25, 4, 31, 1, 22, 3,
        19, 18, 10, 9, 2, 16, 17, 27,
        26, 6, 15, 5, 14, 8, 20, 29,
    )
    assert sorted(m32) == sorted(list(range(1, 32)) + [3])
    expected_m32 = {
        2: ({-2: 341, 0: 771, 2: 128}, Fraction(-213, 620)),
        3: ({-6: 17, -4: 458, -2: 143, 0: 2}, Fraction(-111, 31)),
        4: ({-12: 58, -10: 4}, Fraction(-368, 31)),
    }
    for dimension in (2, 3, 4):
        result = coset_stats(m32, 32, subspaces(32, dimension))
        expected_gains, expected_mean = expected_m32[dimension]
        assert result["best_gains"] == expected_gains
        assert result["mean_best"] == expected_mean
        print("M32", dimension, result)

    d2 = subspaces(64, 2)
    d3 = subspaces(64, 3)
    assert len(d2) == 651
    assert len(d3) == 1395

    random_labels = random_assignment(20260915)
    for dimension, directions in ((2, d2), (3, d3)):
        result = coset_stats(random_labels, 64, directions)
        assert result["positive_cosets"] == (1530 if dimension == 2 else 50)
        assert result["positive_cosets"] == result["positive_triples"]
        assert result["nonuniform_positive"]
        print("M64 random seed=20260915", dimension, result)

    affine, columns, offset = affine_assignment(20260918)
    print("M64 affine columns=", columns, "offset=", offset)
    expected_affine = {
        0: (312, 10), 1: (446, 15), 2: (582, 28),
        4: (728, 21), 8: (1049, 28), 16: (1439, 42),
    }
    for swaps in (0, 1, 2, 4, 8, 16):
        labels = controlled_swaps(affine, 2026091800 + swaps, swaps)
        r2 = coset_stats(labels, 64, d2)
        r3 = coset_stats(labels, 64, d3)
        assert (r2["positive_cosets"], r3["positive_cosets"]) == expected_affine[swaps]
        assert r2["positive_cosets"] == r2["positive_triples"]
        assert r3["positive_cosets"] == r3["positive_triples"]
        assert r2["nonuniform_positive"] and r3["nonuniform_positive"]
        print("M64 affine swaps=", swaps, "d2=", r2, "d3=", r3)

    assert winning_restriction_balance(random_labels, d2) == 60
    assert winning_restriction_balance(random_labels, d3) == 0
    assert winning_restriction_balance(affine, d2) == 17
    assert winning_restriction_balance(affine, d3) == 0
    sixteen = controlled_swaps(affine, 2026091816, 16)
    assert winning_restriction_balance(sixteen, d2) == 73
    assert winning_restriction_balance(sixteen, d3) == 0
    print("exact winning-restriction balance checks: passed")


if __name__ == "__main__":
    main()
