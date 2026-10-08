"""Finite DVR diagnostic for a twelve-row localized-six presentation.

All ranks are exact over F_101.  The 311/30 elementary divisors concern
the selected local-generator presentation at one specified coalition,
not an intrinsic Smith form of the conformal-block kernel lattice.
"""

from collections import Counter
from itertools import product
from math import gcd

import numpy as np

import check_generic_terminal_coherent_syzygies as base
import check_higher_rank_terminal_fusion_grid as fusion


def tensor_support_count(n, d):
    """Sum min(full color blocks, empty color blocks) over all core sets."""
    partial_choices = (1 << d) - 2
    answer = 0
    for states in product((0, 1, 2), repeat=n):
        empty = states.count(0)
        full = states.count(1)
        partial = states.count(2)
        answer += min(empty, full) * partial_choices**partial
    return answer


def check_equal_smith_different_minima():
    """Finite witness for the exact lattice family proved in the note."""
    for n in range(2, 21):
        determinant = n * n + 1
        # Columns (D,0), (0,1) and (D,0), (-n,1) have the same
        # determinant D and first Smith factor gcd(entries)=1.
        assert gcd(determinant, 1) == gcd(determinant, n, 1) == 1
        assert all(
            (x * x + y * y >= determinant)
            for x in range(-n, n + 1)
            for y in range(-n, n + 1)
            if (x, y) != (0, 0) and (x + n * y) % determinant == 0
        )
        assert (-n + n) % determinant == 0


def local_six_slope_ranks():
    prime = base.PRIME
    assert prime == 101
    labels = base.LABELS
    core = set(range(6))
    core_slopes = [i + 1 for i in labels]
    outside_slopes = [i + 17 for i in labels]

    vectors = base.test_vectors()
    determinants = base.triple_determinants(vectors)
    source_eval = np.column_stack([
        np.prod([determinants[triangle] for triangle in column], axis=0)
        % prime for column in base.TRIANGLE_COLUMNS
    ])
    assert base.modular_rank(source_eval) == 462
    assert base.modular_rank(base.twelve_row_pure_core_initial_matrix(6)) == 121

    by_order = [[] for _ in range(4)]
    six_counts = Counter()
    for six in base.SIX_SETS:
        hit = len(core.intersection(six))
        forced_order = max(0, hit - 3)
        six_counts[forced_order] += 1
        gradient = base.gradient_values(six, determinants)
        leading_matchings = []
        for matching in base.MATCHINGS:
            order = 0
            leading = 1
            for left_index, right_index in matching:
                left, right = six[left_index], six[right_index]
                if left in core and right in core:
                    order += 1
                    factor = core_slopes[right] - core_slopes[left]
                elif left in core:
                    factor = outside_slopes[right]
                elif right in core:
                    factor = -outside_slopes[left]
                else:
                    factor = outside_slopes[right] - outside_slopes[left]
                leading = leading * factor % prime
            leading_matchings.append(leading if order == forced_order else 0)
        local_relation = gradient @ np.array(leading_matchings,
                                             dtype=np.int64) % prime
        complement = tuple(label for label in labels if label not in six)
        complement_gradient = base.gradient_values(complement, determinants)
        by_order[forced_order].append(
            local_relation[:, None] * complement_gradient % prime)

    assert tuple(six_counts[i] for i in range(4)) == (662, 225, 36, 1)
    ranks = []
    for order in range(4):
        matrix = np.concatenate(
            [column for group in by_order[:order + 1] for column in group],
            axis=1,
        )
        ranks.append(base.modular_rank(matrix.T))
    assert ranks == [311, 341, 341, 341]
    # The raw columns of positive forced order vanish at t=0; the
    # order-zero columns are exactly the 311-dimensional special image.
    # The normalized order-one columns add 30 independent directions.
    # Hence a maximal minor has valuation exactly 30: rank drop gives
    # the lower bound, and selecting those 30 columns gives the upper.
    return ranks, six_counts


def main():
    check_equal_smith_different_minima()
    ranks, counts = local_six_slope_ranks()
    four_two = fusion.case(4, 2)
    three_four = fusion.case(3, 4)
    assert (four_two["h"], four_two["B"], four_two["threshold"]) == (1, 116, 110)
    assert tensor_support_count(4, 2) == four_two["B"] == 116
    assert tensor_support_count(3, 4) == three_four["threshold"] == 90
    assert (three_four["h"], three_four["B"], three_four["nu"][6]) == (
        341, 171600, 15)
    print("m=12, |T|=6 over F_101: six-set orders", dict(sorted(counts.items())))
    print("normalized local-generator ranks through orders 0,1,2,3:", ranks)
    print("selected presentation Smith slopes: 311 zeros, 30 ones")
    print("fusion numerator order nu_6=15; this is a distinct global scalar")
    print("(4,2): tensor floor=116=B/h; (3,4): tensor floor=90<B/h=171600/341")


if __name__ == "__main__":
    main()
