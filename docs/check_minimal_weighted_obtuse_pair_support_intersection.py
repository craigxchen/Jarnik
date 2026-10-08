"""Exact anchored cut enumeration for pair-support intersection at M<=6.

Every weighted-obtuse sign simplex must pass this combinatorial filter;
the enumeration uses only integer masks, so weight sizes are irrelevant.
"""

from itertools import combinations, permutations


def audit(m):
    r = m - 1
    # Orient each nonconstant unoriented cut to contain row zero.
    cuts = tuple(range(1, 1 << m, 2))[:-1]
    row_pairs = tuple(combinations(range(m), 2))
    distinct = 0
    intersecting = []
    for selected in combinations(cuts, r):
        rows = tuple(sum((1 << j) if (cut >> i) & 1 else 0
                         for j, cut in enumerate(selected)) for i in range(m))
        if len(set(rows)) < m:
            continue
        distinct += 1
        supports = tuple(rows[i] ^ rows[k] for i, k in row_pairs)
        if all(supports[a] & supports[b]
               for a in range(len(supports))
               for b in range(a + 1, len(supports))):
            intersecting.append(selected)
    return distinct, intersecting


if __name__ == "__main__":
    expected = {4: (29, 1), 5: (1015, 0), 6: (126651, 0)}
    for m in expected:
        distinct, good = audit(m)
        assert (distinct, len(good)) == expected[m]
        print("M", m, "distinct-row matrices", distinct,
              "pair-support-intersecting", len(good), "passed")
        if m == 4:
            assert good == [(3, 5, 9)]
    near_miss = (0, 7, 25, 30, 42, 45, 51, 52)
    pair_supports = tuple(a ^ b for a, b in combinations(near_miss, 2))
    assert all(pair_supports[i] & pair_supports[j]
               for i in range(len(pair_supports))
               for j in range(i + 1, len(pair_supports)))
    assert {bin(s).count("1") for s in pair_supports} == {3, 4}
    first_seven = tuple(tuple(1 if not (v >> j) & 1 else -1
                              for j in range(7)) for v in near_miss[:7])
    det = 0
    for perm in permutations(range(7)):
        inversion_count = sum(perm[i] > perm[k]
                              for i in range(7) for k in range(i + 1, 7))
        term = -1 if inversion_count % 2 else 1
        for i in range(7):
            term *= first_seven[i][perm[i]]
        det += term
    assert det == 512
    print("eight-row non-Hadamard support-filter near miss passed")
