"""Exact finite audit of affine-coset characters for one nonlinear Walsh assignment."""

from collections import Counter
from itertools import combinations


LABELS = (
    28, 21, 11, 30, 24, 23, 13, 3,
    12, 7, 25, 4, 31, 1, 22, 3,
    19, 18, 10, 9, 2, 16, 17, 27,
    26, 6, 15, 5, 14, 8, 20, 29,
)
M = 32
B = 5


def parity(n):
    return bin(n).count("1") & 1


def subspaces(dim):
    result = set()
    for basis in combinations(range(1, M), dim):
        elements = {0}
        for vector in basis:
            elements |= {x ^ vector for x in tuple(elements)}
        if len(elements) == 1 << dim:
            result.add(tuple(sorted(elements)))
    return result


def count_restriction_matches(vectors, rows, restriction):
    return sum(
        tuple(parity(LABELS[x] & v) for v in vectors) == restriction
        for x in rows
    )


def exhaustive_counts():
    totals = {}
    for dim in range(1, 6):
        counter = Counter()
        directions = subspaces(dim)
        for vectors in directions:
            restrictions = {
                tuple(parity(a & v) for v in vectors)
                for a in range(1, M)
            }
            restrictions.discard((0,) * len(vectors))
            representatives = {
                min(x ^ v for v in vectors) for x in range(M)
            }
            for x0 in representatives:
                rows = tuple(x0 ^ v for v in vectors)
                for restriction in restrictions:
                    counter[count_restriction_matches(
                        vectors, rows, restriction
                    )] += 1
        totals[dim] = (len(directions), counter)
    return totals


def physical_l1(vectors, x0, ambient_label, b=B):
    """Build all physical coefficients directly, assigning row x to copy 0/1."""
    rows = {x0 ^ v for v in vectors}
    l = {v: parity(ambient_label & v) for v in vectors}
    assert any(l.values())
    c = {
        x0 ^ v: (-1 if l[v] else 1)
        for v in vectors
    }
    copies_used = Counter()
    flipped_copy = {}
    for x, a in enumerate(LABELS):
        flipped_copy[x] = copies_used[a]
        copies_used[a] += 1
    assert max(copies_used.values()) <= b
    total = 0
    for a in range(1, M):
        u = sum(c[x] * (-1 if parity(a & x) else 1)
                for x in rows) // 2
        for copy in range(b):
            v = u
            for x in rows:
                if LABELS[x] == a and flipped_copy[x] == copy:
                    v -= c[x] * (-1 if parity(a & x) else 1)
            total += abs(v)
    return total


def main():
    assert sorted(LABELS) == sorted(tuple(range(1, M)) + (3,))
    expected = {
        1: (31, {0: 103, 1: 274, 2: 119}),
        2: (155, {0: 1001, 1: 1726, 2: 865, 3: 128}),
        3: (155, {0: 1180, 1: 1991, 2: 1020, 3: 147, 4: 2}),
        4: (31, {0: 201, 1: 502, 2: 223, 3: 4}),
        5: (1, {1: 30, 2: 1}),
    }
    totals = exhaustive_counts()
    for dim in range(1, 6):
        n, counts = totals[dim]
        assert (n, dict(counts)) == expected[dim]
        h = 1 << dim
        k = max(counts)
        coefficient_sum = B * M // 2 + h - 2 * k
        unflipped_half = B * (M - 1) / 2
        assert coefficient_sum > unflipped_half
        print("dim=%d h=%d directions=%d max_k=%d min_l1=%d gap=%g" %
              (dim, h, n, k, coefficient_sum,
               coefficient_sum - unflipped_half))

    # The maximal three-match plane is a physical-column fixture.
    plane = (0, 5, 10, 15)
    assert count_restriction_matches(
        plane, plane, tuple(parity(28 & v) for v in plane)
    ) == 3
    assert physical_l1(plane, 0, 28) == B * M // 2 + 4 - 6
    print("physical-column fixture: passed")
    # Translate one nonplanar row to zero; three independent differences
    # make the remaining Walsh signs uniform over these eight patterns.
    for positive in combinations(range(4), 2):
        c = [1 if j in positive else -1 for j in range(4)]
        assert sum(c) == 0
        total = 0
        for pattern in range(8):
            signs = [1] + [(-1) ** ((pattern >> j) & 1) for j in range(3)]
            twice = sum(x * y for x, y in zip(c, signs))
            assert twice % 2 == 0
            total += abs(twice // 2)
        assert total == 6
    print("six balanced nonplanar four-row Fourier identities: passed")


if __name__ == "__main__":
    main()
