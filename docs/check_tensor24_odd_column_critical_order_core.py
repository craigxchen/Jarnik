"""Independent ternary-orientation BFS for the tensor critical-order core."""

import json
from itertools import combinations
from pathlib import Path


ROOT = Path(__file__).resolve().parent
S = json.loads(
    (ROOT / "tensor24_odd_column_certificate_gap_fixture.json").read_text()
)["retained_matrix"]
assert len(S) == 8 and all(len(row) == 22 for row in S)
PATTERN = (1, 1, -1, -1, 1, 1, -1, -1, 1, 1, -1)

certificates = []
for i, j in combinations(range(8), 2):
    c = tuple((S[i][k] - S[j][k]) // 2 for k in range(22))
    if sum(x != 0 for x in c) == 11:
        certificates.append((i, j, c))
assert len(certificates) == 16

supports = tuple(
    sum(1 << j for j, x in enumerate(c) if x)
    for _, _, c in certificates
)


def popcount(value):
    return bin(value).count("1")


def enumerate_states(active):
    """Return exact BFS levels using 0/+1/-1 certificate orientations."""
    ids = tuple(sorted(active))
    index = {certificate: slot for slot, certificate in enumerate(ids)}
    initial = (0,) + (0,) * len(ids)
    levels = [{initial}]
    all_states = {initial}

    for _depth in range(22):
        next_states = set()
        for state in levels[-1]:
            selected = state[0]
            for column in range(22):
                if selected & (1 << column):
                    continue

                required = []
                for certificate in ids:
                    c = certificates[certificate][2]
                    if not (supports[certificate] & (1 << column)):
                        continue
                    position = popcount(selected & supports[certificate])
                    orientation = state[1 + index[certificate]]
                    if position:
                        required.append(
                            c[column] * orientation * PATTERN[position]
                        )
                if required and len(set(required)) != 1:
                    continue

                signs = required[:1] if required else (1, -1)
                for sign in signs:
                    orientations = list(state[1:])
                    for certificate in ids:
                        c = certificates[certificate][2]
                        if (
                            supports[certificate] & (1 << column)
                            and not (selected & supports[certificate])
                        ):
                            orientations[index[certificate]] = c[column] * sign
                    next_states.add(
                        (selected | (1 << column), *orientations)
                    )

        levels.append(next_states)
        all_states.update(next_states)
        if not next_states:
            break

    return levels, all_states


full_levels, full_states = enumerate_states(set(range(16)))
full_counts = [len(level) for level in full_levels]
assert full_counts == [
    1, 44, 214, 200, 326, 408, 264, 300, 372, 396, 216, 120,
    72, 96, 0,
]
assert len(full_states) == 3029
assert len(full_levels) == 15
print("PASS: independent full search has 3,029 states and no 22-column path.")

core_pairs = {
    (1, 4), (2, 5), (2, 6), (2, 7), (3, 5),
    (3, 6), (3, 7), (4, 5), (4, 6), (4, 7),
}
core = {
    certificate
    for certificate, (i, j, _c) in enumerate(certificates)
    if (i, j) in core_pairs
}
assert len(core) == 10
core_levels, core_states = enumerate_states(core)
core_counts = [len(level) for level in core_levels]
assert core_counts == [
    1, 43, 302, 822, 1672, 3202, 5390, 8582, 13098, 18122,
    21838, 23620, 24074, 21628, 15884, 9472, 5466, 3186,
    1124, 200, 28, 0,
]
assert len(core_states) == 177754
assert len(core_levels) == 22
print("PASS: ten-certificate core is infeasible: 177,754 states, no depth-22 path.")
