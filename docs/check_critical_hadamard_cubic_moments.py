"""Exact finite audit of the two-column-deleted Hadamard cubic obstruction.

The proof is in critical_pell_moment_codes.md, Section 8. This checker
does not assert a classification of all length-ten distance-five codes.
"""

from collections import Counter, defaultdict
from fractions import Fraction
from itertools import combinations

from check_one_sided_block_phase_cone import paley


def dot(a, b):
    return sum(x * y for x, y in zip(a, b))


def audit_triple_gaps():
    values = {}
    for a in range(1, 4):
        for b in range(1, 4 - a):
            x = -Fraction(2 * a + b, 3)
            y, z = x + a, x + a + b
            assert x + y + z == 0
            q = x * x + x * y + y * y
            assert q == x * x + x * z + z * z
            assert q == y * y + y * z + z * z
            values[a, b] = q
            if (a, b) == (1, 1):
                assert (x, y, z) == (-1, 0, 1)
            else:
                assert q == Fraction(7, 3)
                assert q not in (1, 4, 9)
    assert set(values) == {(1, 1), (1, 2), (2, 1)}
    for x in range(-7, 8):
        for y in range(-7, 8):
            q, d2 = x * x + x * y + y * y, (x - y) ** 2
            assert 3 * x * y == q - d2
            assert 3 * (x * x + y * y) == 2 * q + d2
    assert 4 ** 2 + 8 ** 2 + 12 ** 2 == 224 > 192


def audit_hadamard():
    H = [[1] + row for row in paley(11)]
    assert all(dot(a, b) == (12 if i == j else 0)
               for i, a in enumerate(H) for j, b in enumerate(H))
    row_sets = signings = deleted_codes = distance_five_codes = 0
    group_types = Counter()
    counts = [bin(mask).count("1") for mask in range(256)]
    for selected in combinations(range(12), 8):
        omitted = [i for i in range(12) if i not in selected]
        orientation = H[omitted[0]]
        T = [[H[i][j] * orientation[j] for j in range(12)] for i in selected]
        K = [[H[i][j] * orientation[j] for j in range(12)] for i in omitted]
        groups = defaultdict(list)
        for j in range(12):
            groups[tuple(row[j] for row in K)].append(j)
        shape = tuple(sorted(map(len, groups.values())))
        assert shape in ((1, 1, 1, 1, 2, 2, 2, 2), (3, 3, 3, 3))
        group_types[shape] += 1
        columns = list(zip(*T))
        pairs = [tuple(group) for group in groups.values() if len(group) == 2]
        differences = [tuple(columns[j][i] - columns[k][i] for i in range(8))
                       for j, k in pairs]
        for (j, k), v in zip(pairs, differences):
            assert dot(v, v) == 24
            assert sum(value != 0 for value in v) == 6
            # Exact row-space projection, with no inverse or floating point.
            assert [sum(v[i] * T[i][c] for i in range(8)) for c in range(12)] == [
                12 * ((c == j) - (c == k)) for c in range(12)]
        assert all(dot(a, b) == 0 for a, b in combinations(differences, 2))
        masks = [sum((entry == 1) << i for i, entry in enumerate(col))
                 for col in columns]
        for signing in range(256):
            rho = [1 if signing >> i & 1 else -1 for i in range(8)]
            u = [dot(rho, v) for v in differences]
            assert all(value in (-12, -8, -4, 0, 4, 8, 12) for value in u)
            assert sum(value * value for value in u) <= 192
            assert set(map(abs, u)) != {0, 4, 8, 12}
            signed = [mask ^ signing for mask in masks]
            defects = [(counts[mask] - 4) ** 2 for mask in signed]
            balanced = [counts[mask] == 4 for mask in signed]
            assert sum(defects) == 24
            for j, k in combinations(range(12), 2):
                a, b = signed[j], signed[k]
                opposite = bool(((a & b) and ((~a & ~b) & 255))
                                or ((a & ~b) and (~a & b)))
                constant = a in (0, 255) or b in (0, 255)
                assert opposite != constant
                if not opposite:
                    defect = 24 - defects[j] - defects[k]
                    count = sum(balanced) - balanced[j] - balanced[k]
                    assert defect <= 8 and defect + count - 20 <= -2
                    distance_five_codes += 1
                deleted_codes += 1
            signings += 1
        row_sets += 1
    assert row_sets == 495 and signings == 126720
    assert deleted_codes == signings * 66
    print(f"PASS: {row_sets} Paley eight-row restrictions, {signings:,} row signings;")
    print(f"{deleted_codes:,} two-column deletions, {distance_five_codes:,} distance-five codes;")
    print(f"pattern groups {dict(group_types)}, exact projections and Bessel bound.")


if __name__ == "__main__":
    audit_triple_gaps()
    audit_hadamard()
    print("PASS: exact cubic-pair invariants, triple gaps, and strict bonus scope.")
