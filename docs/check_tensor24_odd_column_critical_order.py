"""Exact 11-root sign-order obstruction for the fixed 8-by-22 tensor code.

Uses only the literal stored sign matrix and the proved critical odd-moment
root-order lemma. No numerical optimization or coefficient approximation.
"""

from itertools import combinations
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parent
S = json.loads(
    (ROOT / "tensor24_odd_column_certificate_gap_fixture.json").read_text()
)["retained_matrix"]
assert len(S) == 8 and all(len(row) == 22 for row in S)
assert all(entry in (-1, 1) for row in S for entry in row)

# The 11-root lemma forces either PATTERN or its overall negation when
# restricted to the support of a distance-11 row difference.
PATTERN = (1, 1, -1, -1, 1, 1, -1, -1, 1, 1, -1)
certificates = []
for i, j in combinations(range(8), 2):
    c = tuple((S[i][k] - S[j][k]) // 2 for k in range(22))
    if sum(x != 0 for x in c) == 11:
        assert all(x in (-1, 0, 1) for x in c)
        certificates.append((i, j, c))
assert len(certificates) == 16
assert sorted((i, j) for i, j, _ in certificates) == [
    (i, j) for i, j in combinations(range(8), 2)
    if (i in (0, 2, 3, 4)) != (j in (0, 2, 3, 4))
]

supports = tuple(
    sum(1 << j for j, x in enumerate(c) if x) for _, _, c in certificates
)
incident = tuple(
    tuple(k for k, (_, _, c) in enumerate(certificates) if c[j])
    for j in range(22)
)
full_mask = (1 << 22) - 1
seen = set()
depth_counts = [0] * 23


def popcount(mask):
    return bin(mask).count("1")


def search(mask, positive_orientation_bits):
    """Try every next column and sign, merging equivalent exact prefixes.

    An untouched certificate has no orientation yet. Once touched, its
    orientation is stored as one bit. Its pattern position is exactly the
    number of already selected columns in its support, so the selected mask
    and orientation bits contain all information needed for the suffix.
    """
    if mask == full_mask:
        return ()
    state = (mask, positive_orientation_bits)
    if state in seen:
        return None
    seen.add(state)
    depth_counts[popcount(mask)] += 1

    for j in range(22):
        if mask & (1 << j):
            continue
        required_sign = 0
        for k in incident[j]:
            position = popcount(mask & supports[k])
            if position == 0:
                continue
            orientation = 1 if positive_orientation_bits & (1 << k) else -1
            sign = orientation * PATTERN[position] * certificates[k][2][j]
            if required_sign and required_sign != sign:
                break
            required_sign = sign
        else:
            possibilities = (required_sign,) if required_sign else (1, -1)
            for sign in possibilities:
                new_bits = positive_orientation_bits
                for k in incident[j]:
                    if mask & supports[k] == 0:
                        if certificates[k][2][j] * sign == 1:
                            new_bits |= 1 << k
                suffix = search(mask | (1 << j), new_bits)
                if suffix is not None:
                    return ((j, sign),) + suffix
    return None


answer = search(0, 0)
assert answer is None, answer
assert len(seen) == 3029, len(seen)
assert sum(depth_counts) == len(seen)
print("PASS: 16 distance-11 certificates; no common order and signs")
print("Exact memoized states:", len(seen))
print("States by selected-column count:", depth_counts)
