#!/usr/bin/env python3
"""Full-rank certificate on the actual q_i^2/D_i balanced-cut boundary."""

from itertools import combinations
import random

import numpy as np


PRIME = 1009
RANDOM = random.Random(20260913)
CUTS = [set(cut) for cut in combinations(range(8), 4) if 0 in cut]
TRIPLES = list(combinations(range(7), 3))
QUADRATIC_PAIRS = [(i, j) for i in range(35) for j in range(i, 35)]


def inverse(value):
    return pow(int(value) % PRIME, PRIME - 2, PRIME)


def det3(a, b, c):
    return (
        a[0] * (b[1] * c[2] - b[2] * c[1])
        - a[1] * (b[0] * c[2] - b[2] * c[0])
        + a[2] * (b[0] * c[1] - b[1] * c[0])
    ) % PRIME


def boundary_evaluation():
    while True:
        cut = RANDOM.choice(CUTS)
        complement = set(range(8)) - cut
        x = dict(zip(sorted(cut), RANDOM.sample(range(1, PRIME), 4)))
        y = dict(zip(sorted(complement), RANDOM.sample(range(1, PRIME), 4)))
        coefficient_a, coefficient_b, coefficient_c = [
            RANDOM.randrange(1, PRIME) for _ in range(3)
        ]

        scaling = {}
        for i in cut:
            same_side_difference = 1
            for k in cut:
                if k != i:
                    same_side_difference *= x[k] - x[i]
            scaling[i] = (
                (coefficient_a * x[i] + coefficient_b) ** 2
                * inverse(x[i] * same_side_difference)
            ) % PRIME

        x_product = 1
        for value in x.values():
            x_product = x_product * value % PRIME
        for j in complement:
            reciprocal_difference = 1
            for k in complement:
                if k != j:
                    reciprocal_difference = (
                        reciprocal_difference * (inverse(y[k]) - inverse(y[j]))
                    ) % PRIME
            scaling[j] = (
                (coefficient_c + coefficient_b * inverse(y[j])) ** 2
                * inverse(y[j] * x_product * reciprocal_difference)
            ) % PRIME

        if all(scaling.values()):
            break

    rows = [
        (scaling[i], scaling[i] * x[i] % PRIME, 0)
        if i in cut
        else (scaling[i], 0, scaling[i] * y[i] % PRIME)
        for i in range(8)
    ]
    assert all(sum(row[column] for row in rows) % PRIME == 0 for column in range(3))
    pluckers = np.array(
        [det3(rows[i], rows[j], rows[k]) for i, j, k in TRIPLES],
        dtype=np.int64,
    )
    return np.array(
        [int(pluckers[i]) * int(pluckers[j]) % PRIME for i, j in QUADRATIC_PAIRS],
        dtype=np.int64,
    )


class EchelonBasis:
    def __init__(self):
        self.rows = {}

    def add(self, row):
        row = row.copy()
        for pivot, basis_row in self.rows.items():
            if row[pivot]:
                row = (row - int(row[pivot]) * basis_row) % PRIME
        nonzero = np.flatnonzero(row)
        if len(nonzero) == 0:
            return False
        pivot = int(nonzero[0])
        row = row * inverse(row[pivot]) % PRIME
        self.rows[pivot] = row
        self.rows = dict(sorted(self.rows.items()))
        return True


def hook_dimension_222(dimension):
    partition = (2, 2, 2)
    numerator = 1
    denominator = 1
    for i, row_length in enumerate(partition):
        for j in range(row_length):
            hook = row_length - j + sum(length > j for length in partition[i + 1 :])
            numerator *= dimension + j - i
            denominator *= hook
    assert numerator % denominator == 0
    return numerator // denominator


def main():
    target = hook_dimension_222(7)
    assert target == 490
    basis = EchelonBasis()
    for sample in range(target):
        assert basis.add(boundary_evaluation()), f"rank stalled at sample {sample}"
    assert len(basis.rows) == target
    print("PASS: 490 actual balanced-boundary evaluations have rank 490 over F_1009.")


if __name__ == "__main__":
    main()
