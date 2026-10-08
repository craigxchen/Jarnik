#!/usr/bin/env python3
"""Exact F2 certificate for H^1(S_12, even subsets / complements)."""

from __future__ import annotations

from collections import deque
from functools import reduce
from hashlib import sha256
from typing import Optional


N_COORD = 10
N_GENERATORS = 11
N_VARIABLES = N_COORD*N_GENERATORS
ALL_TWELVE = (1 << 12)-1
IDENTITY = [1 << j for j in range(N_COORD)]
EXPECTED_FREE = [1, 12, 23, 34, 45, 56, 67, 78, 89, 99, 109]
EXPECTED_RREF_SHA256 = "f49378980e1deaebda82fd425ad983aff39f3aec87b253208e0f8e274f8ca68e"


def parity(x: int) -> int:
    return bin(x).count("1") & 1


def lift(x: int) -> int:
    """Lift ten coordinates to the even 12-bit representative with bit 11 zero."""
    return x | (parity(x) << 10)


def generator_action(i: int, x: int) -> int:
    """Apply the transposition (i,i+1), then choose bit 11 zero."""
    y = lift(x)
    if ((y >> i) & 1) != ((y >> (i+1)) & 1):
        y ^= (1 << i) | (1 << (i+1))
    if (y >> 11) & 1:
        y ^= ALL_TWELVE
    return y & ((1 << N_COORD)-1)


ACTION = [
    [generator_action(i, 1 << j) for j in range(N_COORD)]
    for i in range(N_GENERATORS)
]


def mat_vec(matrix: list[int], x: int) -> int:
    answer = 0
    for j, column in enumerate(matrix):
        if (x >> j) & 1:
            answer ^= column
    return answer


def mat_mul(left: list[int], right: list[int]) -> list[int]:
    return [mat_vec(left, column) for column in right]


def mat_sum(*matrices: list[int]) -> list[int]:
    return [reduce(int.__xor__, columns, 0) for columns in zip(*matrices)]


def add_vector_equation(rows: list[int], terms: list[tuple[int, list[int]]]) -> None:
    """Append the ten scalar rows of sum matrix*X_block=0."""
    for output in range(N_COORD):
        row = 0
        for block, matrix in terms:
            for input_coordinate, column in enumerate(matrix):
                if (column >> output) & 1:
                    row ^= 1 << (N_COORD*block+input_coordinate)
        rows.append(row)


def relation_matrix(number_of_generators: int = N_GENERATORS) -> list[int]:
    rows: list[int] = []
    for i in range(number_of_generators):
        add_vector_equation(rows, [(i, mat_sum(IDENTITY, ACTION[i]))])
    for i in range(number_of_generators):
        for j in range(i+2, number_of_generators):
            add_vector_equation(rows, [
                (i, mat_sum(IDENTITY, ACTION[j])),
                (j, mat_sum(IDENTITY, ACTION[i])),
            ])
    for i in range(number_of_generators-1):
        j = i+1
        add_vector_equation(rows, [
            (i, mat_sum(IDENTITY, mat_mul(ACTION[i], ACTION[j]), ACTION[j])),
            (j, mat_sum(ACTION[i], IDENTITY, mat_mul(ACTION[j], ACTION[i]))),
        ])
    return rows


def rref(rows: list[int], columns: int) -> tuple[list[int], list[int]]:
    work = [row for row in rows if row]
    pivots: list[int] = []
    rank = 0
    for column in range(columns):
        pivot = next((r for r in range(rank, len(work))
                      if (work[r] >> column) & 1), None)
        if pivot is None:
            continue
        work[rank], work[pivot] = work[pivot], work[rank]
        for r in range(len(work)):
            if r != rank and ((work[r] >> column) & 1):
                work[r] ^= work[rank]
        pivots.append(column)
        rank += 1
    return work[:rank], pivots


def satisfies(rows: list[int], vector: int) -> bool:
    return all(parity(row & vector) == 0 for row in rows)


def coboundary(coordinate: int) -> int:
    answer = 0
    b = 1 << coordinate
    for i in range(N_GENERATORS):
        answer |= (generator_action(i, b) ^ b) << (N_COORD*i)
    return answer


def restricted_coboundary(coordinate: int, number_of_generators: int) -> int:
    answer = 0
    b = 1 << coordinate
    for i in range(number_of_generators):
        answer |= (generator_action(i, b) ^ b) << (N_COORD*i)
    return answer


def orbit_sizes(translations: list[int]) -> list[int]:
    unseen = set(range(1 << N_COORD))
    sizes: list[int] = []
    while unseen:
        start = min(unseen)
        orbit = {start}
        queue = deque([start])
        while queue:
            x = queue.popleft()
            for i in range(N_GENERATORS):
                y = generator_action(i, x) ^ translations[i]
                if y not in orbit:
                    orbit.add(y)
                    queue.append(y)
        unseen.difference_update(orbit)
        sizes.append(len(orbit))
    return sorted(sizes)


def permutation_action(permutation: list[int], x: int) -> int:
    y = lift(x)
    image = 0
    for i in range(12):
        if (y >> i) & 1:
            image |= 1 << permutation[i]
    if (image >> 11) & 1:
        image ^= ALL_TWELVE
    return image & ((1 << N_COORD)-1)


def linear_orbit_sizes(actions: list[list[int]]) -> list[int]:
    unseen = set(range(1 << N_COORD))
    sizes: list[int] = []
    while unseen:
        start = min(unseen)
        orbit = {start}
        queue = deque([start])
        while queue:
            x = queue.popleft()
            for action in actions:
                y = mat_vec(action, x)
                if y not in orbit:
                    orbit.add(y)
                    queue.append(y)
        unseen.difference_update(orbit)
        sizes.append(len(orbit))
    return sorted(sizes)


def partitions(total: int, maximum: Optional[int] = None) -> list[tuple[int, ...]]:
    if total == 0:
        return [()]
    if maximum is None or maximum > total:
        maximum = total
    answer: list[tuple[int, ...]] = []
    for first in range(maximum, 0, -1):
        for rest in partitions(total-first, first):
            answer.append((first,)+rest)
    return answer


def main() -> None:
    # Verify the type-A_11 Coxeter relations on the module itself.
    for i in range(N_GENERATORS):
        assert mat_mul(ACTION[i], ACTION[i]) == IDENTITY
    for i in range(N_GENERATORS):
        for j in range(i+2, N_GENERATORS):
            assert mat_mul(ACTION[i], ACTION[j]) == mat_mul(ACTION[j], ACTION[i])
    for i in range(N_GENERATORS-1):
        left = mat_mul(mat_mul(ACTION[i], ACTION[i+1]), ACTION[i])
        right = mat_mul(mat_mul(ACTION[i+1], ACTION[i]), ACTION[i+1])
        assert left == right

    rows = relation_matrix()
    assert len(rows) == 660
    reduced, pivots = rref(rows, N_VARIABLES)
    free = [column for column in range(N_VARIABLES) if column not in set(pivots)]
    assert len(pivots) == 99
    assert free == EXPECTED_FREE
    digest = sha256(b"".join(row.to_bytes(14, "little") for row in reduced)).hexdigest()
    assert digest == EXPECTED_RREF_SHA256

    boundaries = [coboundary(j) for j in range(N_COORD)]
    assert len(rref(boundaries, N_VARIABLES)[1]) == 10
    assert all(satisfies(rows, boundary) for boundary in boundaries)

    # Class of the odd-subset torsor based at {11}.  In canonical coordinates
    # [{10,11}] is represented by its complement {0,...,9}.
    odd = ((1 << N_COORD)-1) << (N_COORD*10)
    assert satisfies(rows, odd)
    assert len(rref(boundaries+[odd], N_VARIABLES)[1]) == 11
    # Since dim Z^1=11, these ten boundaries and odd span every cocycle.

    zero_translations = [0]*N_GENERATORS
    odd_translations = [0]*N_GENERATORS
    odd_translations[10] = (1 << N_COORD)-1
    assert orbit_sizes(zero_translations) == [1, 66, 462, 495]
    assert orbit_sizes(odd_translations) == [12, 220, 792]

    # The stabilizer of label 11 is S_11, generated by s_0,...,s_9.
    # On V this is the even-subset module on the remaining eleven labels.
    s11_rows = relation_matrix(10)
    assert len(s11_rows) == 550
    assert len(rref(s11_rows, 100)[1]) == 90
    s11_boundaries = [restricted_coboundary(j, 10) for j in range(N_COORD)]
    assert len(rref(s11_boundaries, 100)[1]) == 10
    assert all(satisfies(s11_rows, boundary) for boundary in s11_boundaries)
    assert orbit_sizes([0]*N_GENERATORS) == [1, 66, 462, 495]
    unseen = set(range(1 << N_COORD))
    s11_sizes: list[int] = []
    while unseen:
        start = min(unseen)
        orbit = {start}
        queue = deque([start])
        while queue:
            x = queue.popleft()
            for i in range(10):
                y = generator_action(i, x)
                if y not in orbit:
                    orbit.add(y)
                    queue.append(y)
        unseen.difference_update(orbit)
        s11_sizes.append(len(orbit))
    assert sorted(s11_sizes) == [1, 11, 55, 165, 330, 462]

    # Point stabilizers for subset weights k=0,...,6.  For k=6 include
    # the permutation interchanging the two six-element blocks.
    expected_counts = [4, 6, 9, 10, 12, 12]
    for k, expected in enumerate(expected_counts):
        generators = [ACTION[i] for i in range(N_GENERATORS)
                      if k == 0 or i != k-1]
        assert len(linear_orbit_sizes(generators)) == expected
    half_swap = list(range(12))
    for i in range(6):
        half_swap[i], half_swap[i+6] = i+6, i
    half_swap_matrix = [permutation_action(half_swap, 1 << j)
                        for j in range(N_COORD)]
    generators = [ACTION[i] for i in range(N_GENERATORS) if i != 5]
    generators.append(half_swap_matrix)
    assert len(linear_orbit_sizes(generators)) == 10

    # Every simultaneous stabilizer of a target subset class and a V-class
    # contains an odd transposition: their four membership cells total 12.
    # Thus intersecting a target stabilizer with A_12 does not split an orbit.
    for k in range(7):
        target = set(range(k))
        for x in range(1 << N_COORD):
            subset = {i for i in range(12) if (lift(x) >> i) & 1}
            cells = [
                {i for i in range(12) if (i in target) == target_side
                 and (i in subset) == subset_side}
                for target_side in (False, True)
                for subset_side in (False, True)
            ]
            assert max(map(len, cells)) >= 3

    # Any permutation of twelve letters has at least six cycles on pairs.
    from math import gcd
    pair_cycle_counts = []
    for cycle_lengths in partitions(12):
        count = sum(length//2 for length in cycle_lengths)
        count += sum(gcd(cycle_lengths[i], cycle_lengths[j])
                     for i in range(len(cycle_lengths))
                     for j in range(i+1, len(cycle_lengths)))
        pair_cycle_counts.append(count)
    assert min(pair_cycle_counts) == 6

    print("PASS: 660 Coxeter cocycle rows have rank 99, so dim Z1=11.")
    print("PASS: dim B1=10 and the odd-subset cocycle generates H1 of dimension one.")
    print("PASS: linear orbit sizes 1,66,495,462 and odd-torsor sizes 12,220,792.")
    print("PASS: for the S11 point stabilizer, H1=0 and there are six linear orbits.")
    print("PASS: subset-stabilizer orbit counts are 4,6,9,10,12,12,10; pair-cycle minimum is 6.")


if __name__ == "__main__":
    main()
