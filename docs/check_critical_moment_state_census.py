"""Exact compressed-state census by linear extensions, without graph traversal."""
from collections import Counter
from itertools import product
from random import Random

from check_critical_moment_compressed_automaton import Automaton, PAIRS


def predecessors(machine, heights):
    pred = [0] * 8
    signs = [(-1) ** (machine.delta[i] - heights[i]) for i in range(8)]
    for (i, j), distance, kind, orientation, coefficient_sign in zip(
        PAIRS, machine.distance, machine.kind, machine.orientation, machine.kappa_sign
    ):
        height = orientation * (heights[i] - heights[j])
        allowed = []
        for mismatch in (0, 1):
            phase = (distance + heights[i] + heights[j]
                     - machine.delta[i] - machine.delta[j] + 2 * mismatch) % 4
            valid = (
                (kind == "C" and height == (0, 1, 2, 1)[phase])
                or (kind == "A" and height == (0, 1, 0, -1)[phase])
                or (kind == "B" and ((height == 0 and phase == 0)
                                     or height == (2, 1, 2, 3)[phase]))
            )
            comparison = coefficient_sign * (-1) ** mismatch
            if signs[i] != signs[j]:
                valid &= comparison == (1 if signs[i] > signs[j] else -1)
            if valid:
                allowed.append(comparison)
        if not allowed:
            return None
        if len(allowed) == 1:
            lower, upper = (j, i) if allowed[0] == 1 else (i, j)
            pred[upper] |= 1 << lower
    return tuple(pred)


def extension_count(pred):
    if pred is None:
        return 0
    dp = [0] * 256
    dp[0] = 1
    for selected in range(256):
        for row in range(8):
            if not selected >> row & 1 and not pred[row] & ~selected:
                dp[selected | 1 << row] += dp[selected]
    return dp[-1]


def census(machine, audit=False):
    a_heights = ({(0,) + a for a in product((0, 1), repeat=3)}
                 | {(0,) + a for a in product((0, -1), repeat=3)})
    assert len(a_heights) == 15
    histogram = Counter()
    rng = Random(202609210 + machine.q + 5 * machine.s_parity)
    for a in sorted(a_heights):
        for b in product(*[(0, d, 2 * d) for d in machine.delta[4:]]):
            heights = a + b
            pred = predecessors(machine, heights)
            count = extension_count(pred)
            if count:
                histogram[count] += 1
            if audit:
                # Independent calls to the original state validator check the
                # conversion of pair constraints to predecessor relations.
                for _ in range(8):
                    order = list(range(8))
                    rng.shuffle(order)
                    selected = 0
                    accepted = pred is not None
                    if accepted:
                        for row in order:
                            if pred[row] & ~selected:
                                accepted = False
                                break
                            selected |= 1 << row
                    assert accepted == machine.valid((heights, tuple(order)))
    return histogram


if __name__ == "__main__":
    base = Counter({576: 70, 720: 112, 1440: 56, 5040: 16, 40320: 1})
    expected = [base, base + Counter({144: 1, 240: 3, 720: 3}),
                base + Counter({96: 1, 240: 4, 720: 4})]
    expected += [expected[1], expected[0]]
    for parity in (0, 1):
        for q in range(5):
            histogram = census(Automaton(parity, q), audit=True)
            assert histogram == expected[q]
            total = sum(n * multiplicity for n, multiplicity in histogram.items())
            print(f"PASS: parity={parity}, q={q}: {sum(histogram.values())} "
                  f"height vectors, {total} valid states.")
    print("PASS: exact subset-DP counts and 97,200 reference-validator comparisons.")
