#!/usr/bin/env python3
"""Exact checks for semibent_eight_point_phase_obstruction.md.

The proof is in the note. This checks all matrix identities, every small
nonunit-support pattern used in the counting step, and integer examples
of the elementary phase bounds and the twisted Pell obstruction.
"""

from itertools import combinations
from math import gcd


def bits(x):
    return tuple((x >> j) & 1 for j in range(3))


labels = list(range(8))
H = []
for x in labels:
    xx = bits(x)
    H.append([
        (-1) ** (xx[0] * xx[1] + sum(a * b for a, b in zip(xx, bits(y))))
        for y in labels
    ])
A = [[(1 - h) // 2 for h in row] for row in H]
assert A[0] == [0] * 8
assert all(sum(H[x][y] * H[z][y] for y in labels) == 8 * (x == z)
           for x in labels for z in labels)
S = [sum(H[x][y] for x in labels) for y in labels]
U = [y for y in labels if S[y] != 0]
B = [y for y in labels if S[y] == 0]
assert len(U) == len(B) == 4
assert sorted(S[y] for y in U) == [-4, 4, 4, 4]

for y in B:
    actual = [sum(-H[x][y] * A[x][j] for x in range(1, 8)) for j in labels]
    assert actual == [4 * (j == y) for j in labels]

for y, z in combinations(U, 2):
    epsilon = S[y] // S[z]
    coefficients = [-(H[x][y] - epsilon * H[x][z]) // 2 for x in range(1, 8)]
    assert all(c in (-1, 0, 1) for c in coefficients)
    assert sum(abs(c) for c in coefficients) <= 4
    actual = [sum(coefficients[x - 1] * A[x][j] for x in range(1, 8))
              for j in labels]
    assert actual == [2 * ((j == y) - epsilon * (j == z)) for j in labels]

support_checks = 0
for nb in range(3):
    for balanced in combinations(B, nb):
        for nu in range(2):
            for unbalanced in combinations(U, nu):
                selected = balanced + unbalanced
                patterns = {tuple(A[x][j] for j in selected) for x in labels}
                assert len(patterns) <= 6
                support_checks += 1
assert support_checks == 55

phase_checks = 0
for a in range(1, 101):
    for b in range(1, 101):
        if gcd(a, b) != 1 or (a + b) % 2 == 0:
            continue
        norm = a * a + b * b
        assert (2 * a * b) ** 2 >= 2 * norm
        assert (4 * a * b * (a * a - b * b)) ** 2 >= norm ** 3
        phase_checks += 1

x, b = 1, 0
for _ in range(12):
    x, b = 9 * x + 20 * b, 4 * x + 9 * b
    a = x - 2 * b
    assert x * x - 5 * b * b == 1
    assert a * a + 4 * a * b - b * b == 1
    assert gcd(a, b) == 1 and (a + b) % 2 == 1

assert 117649 < 29 * 4096
assert 29 < 392

spectrum_counts = {}
general_identity_checks = 0
for switch in range(128):
    matrix = []
    for x in labels:
        sign = 1 if x == 0 else (-1) ** ((switch >> (x - 1)) & 1)
        matrix.append([sign * (-1) ** sum(a * b for a, b in zip(bits(x), bits(y)))
                       for y in labels])
    incidence = [[(1 - h) // 2 for h in row] for row in matrix]
    sums = [sum(matrix[x][y] for x in labels) for y in labels]
    assert sum(s * s for s in sums) == 64
    assert len({s % 4 for s in sums}) == 1
    spectrum = tuple(sorted(map(abs, sums), reverse=True))
    spectrum_counts[spectrum] = spectrum_counts.get(spectrum, 0) + 1
    for y in labels:
        if sums[y] != 0:
            continue
        actual = [sum(-matrix[x][y] * incidence[x][j] for x in range(1, 8))
                  for j in labels]
        assert actual == [4 * (j == y) for j in labels]
        general_identity_checks += 1
    for y, z in combinations(labels, 2):
        if sums[y] == 0 or abs(sums[y]) != abs(sums[z]):
            continue
        epsilon = sums[y] // sums[z]
        coefficients = [-(matrix[x][y] - epsilon * matrix[x][z]) // 2
                        for x in range(1, 8)]
        assert sum(abs(c) for c in coefficients) <= 4
        actual = [sum(coefficients[x - 1] * incidence[x][j] for x in range(1, 8))
                  for j in labels]
        assert actual == [2 * ((j == y) - epsilon * (j == z)) for j in labels]
        general_identity_checks += 1
assert spectrum_counts == {
    (8, 0, 0, 0, 0, 0, 0, 0): 8,
    (4, 4, 4, 4, 0, 0, 0, 0): 56,
    (6, 2, 2, 2, 2, 2, 2, 2): 64,
}
assert general_identity_checks == 1960
print("Passed: 10 exact phase-projection identities; 55 support patterns; "
      f"{phase_checks} primitive Gaussian phase checks; 12 twisted Pell identities.")
print("Passed: all 128 Walsh row switchings; spectra counts 8, 56, 64; "
      "1,960 additional exact projection identities.")
