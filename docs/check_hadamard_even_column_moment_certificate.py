"""Exact eight-row integer-moment certificates and the complete Paley24 audit.

The infinite statements are proved in critical_pell_moment_codes.md, Section9.
These fixtures are sign codes; no Gaussian endpoint realization is asserted.
"""

import json
from collections import Counter
from itertools import combinations, product
from pathlib import Path

from check_one_sided_block_phase_cone import paley


def dot(a, b):
    return sum(x * y for x, y in zip(a, b))


def retained_certificate(T, deleted, chosen):
    L = len(T[0])
    columns = list(zip(*T))
    assert all(sum(x == 1 for x in col) % 2 == 0 for col in columns)
    lam = columns[chosen]
    assert sum(lam) == 0 and chosen not in deleted
    retained = [j for j in range(L) if j not in deleted]
    numerators = [dot(lam, columns[j]) for j in retained]
    assert all(n % 4 == 0 for n in numerators)
    c = [n // 4 for n in numerators]
    assert c[retained.index(chosen)] == 2
    assert sum(x * x for x in c) <= L // 2
    assert 0 < sum(map(abs, c)) <= sum(x * x for x in c) - 2 <= L // 2 - 2
    return c


def parity_padded_audit(S):
    m = len(S[0])
    s = (m + 1) // 4
    padding, L = 4 * s + 3 - m, 4 * s + 4
    assert 1 <= padding <= 4
    assert min(sum(a != b for a, b in zip(x, y))
               for x, y in combinations(S, 2)) >= 2 * s + 1
    lam = []
    for row in S:
        p = 1
        for value in row:
            p *= value
        lam.append(p)
    assert sum(lam) == 0
    assert all(sum(row[j] == 1 for row in S) % 2 == 0 for j in range(m))
    T = [row + [parity] + [1] * padding for row, parity in zip(S, lam)]
    A = [[0] * 8 for _ in range(8)]
    for i, j in combinations(range(8), 2):
        value = dot(T[i], T[j])
        assert value <= 0 and value % 4 == 0
        A[i][j] = A[j][i] = -value // 4
    K = sum(A[i][j] for i, j in combinations(range(8), 2))
    form = sum(lam[i] * A[i][j] * lam[j] for i in range(8) for j in range(8))
    nums = [sum(lam[i] * S[i][j] for i in range(8)) for j in range(m)]
    assert all(n % 4 == 0 for n in nums)
    c = [n // 4 for n in nums]
    energy = sum(x * x for x in c)
    assert 4 * energy == 2 * L - 16 - form
    assert abs(form) <= 2 * K
    assert L - 8 - K <= 2 * energy <= L - 8 + K
    if K <= 4 and L - K > 8:
        assert 0 < sum(map(abs, c)) <= 2 * s
    return K


def audit_weighted_newton_small():
    # Nonzero distinct absolute nodes; multiplicity two is included explicitly.
    nodes = (1, -2, 3, -4)
    for c in product(range(-2, 3), repeat=4):
        mass = sum(map(abs, c))
        if not mass:
            continue
        for s in (1, 2, 3, 4):
            if mass <= 2 * s:
                assert any(sum(v * a ** (2 * r - 1) for v, a in zip(c, nodes))
                           for r in range(1, s + 1))


def main():
    fixture = json.loads(Path(__file__).with_name(
        "paley24_even_column_moment_fixture.json").read_text())
    H = [[1] + row for row in paley(23)]
    assert H == fixture["hadamard"]
    assert all(dot(a, b) == (24 if i == j else 0)
               for i, a in enumerate(H) for j, b in enumerate(H))
    column_masks = [sum((H[i][j] > 0) << i for i in range(24)) for j in range(24)]
    counts = [bin(i).count("1") for i in range(4096)]
    histogram, positive = Counter(), []
    selected_count = code_count = 0
    for rows in combinations(range(24), 8):
        mask = sum(1 << i for i in rows)
        plus_counts = []
        for col in column_masks:
            value = col & mask
            plus_counts.append(counts[value & 4095] + counts[value >> 12])
        balanced = plus_counts.count(4)
        histogram[balanced] += 1
        selected_count += 1
        if balanced < 14:
            continue
        assert sum(r in (0, 8) for r in plus_counts) == 1
        positive_deletions = [j for j in range(1, 24)
                              if balanced - 12 - (plus_counts[j] - 4) ** 2
                              - (plus_counts[j] == 4) > 0]
        if not positive_deletions:
            continue
        assert balanced == 15
        assert Counter(min(r, 8 - r) for r in plus_counts) == {0: 1, 2: 8, 4: 15}
        assert positive_deletions == [j for j, r in enumerate(plus_counts) if r == 4]
        T = [H[i] for i in rows]
        columns = list(zip(*T))
        for deleted in positive_deletions:
            keep = [j for j in range(24) if j not in (0, deleted)]
            lam = columns[deleted]
            assert sum(lam) == 0
            nums = [dot(lam, columns[j]) for j in keep]
            assert all(v % 4 == 0 for v in nums)
            c = [v // 4 for v in nums]
            assert set(c) <= {-1, 0, 1} and sum(map(abs, c)) == 8
            assert sum(v * v for v in c) == 8
            chosen = next(j for j in positive_deletions if j != deleted)
            retained_certificate(T, {0, deleted}, chosen)
            S = [[row[j] for j in keep] for row in T]
            assert min(sum(a != b for a, b in zip(x, y))
                       for x, y in combinations(S, 2)) == 11
            assert parity_padded_audit(S) == 0
            code_count += 1
        positive.append(list(rows))
    assert selected_count == fixture["enumerated_row_sets"] == 735471
    assert dict(histogram) == {int(k): v for k, v in fixture["balanced_count_histogram"].items()}
    assert positive == fixture["positive_row_sets"] and len(positive) == 759
    assert code_count == fixture["positive_deleted_codes"] == 11385
    example = fixture["representative"]
    assert len(set(example["deleted_columns"])) == 2
    assert example["retained_columns"] == [
        j for j in range(24) if j not in example["deleted_columns"]]
    assert example["rows"] in positive
    T = [H[i] for i in example["rows"]]
    lam = example["zero_sum_row_coefficients"]
    assert sum(lam) == 0
    assert [sum(lam[i] * T[i][j] for i in range(8)) for j in example["retained_columns"]] == [
        4 * c for c in example["column_coefficients"]]

    # A literal parity-padded code with positive Gram defect, outside K<=4.
    # Its identity is still exact; the sufficient small-defect bound is not used.
    labels = list(range(1, 8)) * 2 + [1, 2, 4, 3]
    S = [[1 if bin(x & a).count("1") % 2 == 0 else -1 for a in labels]
         for x in range(8)]
    assert parity_padded_audit(S) == 12
    audit_weighted_newton_small()
    print(f"PASS: all {selected_count:,} Paley24 row sets; {len(positive)} positive profiles,")
    print(f"{code_count:,} exact deleted- and retained-column moment certificates;")
    print("parity-padded energy identities at K=0 and K=12, weighted Newton fixtures.")


if __name__ == "__main__":
    main()
