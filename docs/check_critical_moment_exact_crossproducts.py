"""Exact 31-block flow audit for the first-moment/cross-product relaxation.

This checker does not materialize the giant Euler word or solve the odd
moments beyond degree one. The continuous block-shift step is proved in the
companion note.
"""

from collections import Counter
from fractions import Fraction
import json
from pathlib import Path

from check_critical_moment_compressed_automaton import Automaton, PAIRS


def solve_many(matrix, targets):
    """Solve a square rational system against several column targets."""
    n = len(matrix)
    q = len(targets)
    rows = [
        [Fraction(x) for x in matrix[i]]
        + [Fraction(targets[j][i]) for j in range(q)]
        for i in range(n)
    ]
    for col in range(n):
        pivot = next(row for row in range(col, n) if rows[row][col])
        rows[col], rows[pivot] = rows[pivot], rows[col]
        scale = rows[col][col]
        rows[col] = [x / scale for x in rows[col]]
        for row in range(n):
            if row == col:
                continue
            factor = rows[row][col]
            if factor:
                rows[row] = [x - factor * y for x, y in zip(rows[row], rows[col])]
    return [[rows[i][n + j] for i in range(n)] for j in range(q)]


def run():
    data = json.loads((Path(__file__).with_name(
        "critical_moment_compressed_graph_even_q2_summary.json")).read_text())
    machine = Automaton(0, 2)
    walks = data["active_cycle_witnesses"]
    correction = data["integer_signed_cycle_correction"]["coefficients"]
    circ = data["active_only_integer_flow_per_edge_by_cut_size"]
    assert len(walks) == len(correction) == len(PAIRS) + 1 == 29
    assert circ == [0, 10, 28, 56, 70, 56, 28, 10, 0]

    cross = [(i, k) for i in range(4) for k in range(4, 8)]
    ref = (0, 4)
    other = [pair for pair in cross if pair != ref]
    assert len(other) == 15
    label_rows = [[walk["labels"][row] for walk in walks] for row in range(29)]
    targets = []
    for pair in other:
        target = [0] * 29
        target[1 + PAIRS.index(pair)] = 4
        targets.append(target)
    basis = solve_many(label_rows, targets)
    assert all(x.denominator == 1 for column in basis for x in column)
    basis = [[int(x) for x in column] for column in basis]
    assert max(abs(x) for column in basis for x in column) <= 6
    for j, column in enumerate(basis):
        expected = targets[j]
        actual = [sum(label_rows[row][k] * column[k] for k in range(29))
                  for row in range(29)]
        assert actual == expected

    walk_edges = []
    for walk in walks:
        state = machine.terminal
        edges = []
        derived = [0] * 29
        for mask in walk["path_masks"]:
            assert mask not in (0, 255)
            transitions = {m: nxt for m, nxt, _ in machine.transitions(state)}
            assert mask in transitions
            edges.append((state, mask))
            derived[0] += 1
            for r, (i, k) in enumerate(PAIRS, 1):
                derived[r] += ((mask >> i) & 1) != ((mask >> k) & 1)
            state = transitions[mask]
        assert state == machine.terminal
        assert derived == walk["labels"]
        walk_edges.append(Counter(edges))

    def signed_flow(coeffs):
        out = Counter()
        for coefficient, edges in zip(coeffs, walk_edges):
            for edge, count in edges.items():
                out[edge] += coefficient * count
        return out

    # e_l is the cumulative cross-label vector after closed block l.
    prefixes = [[0] * 15]
    for j in range(15):
        plus = [0] * 15
        minus = [0] * 15
        plus[j], minus[j] = 4, -4
        prefixes.extend((plus, minus))
    prefixes.append([0] * 15)
    assert len(prefixes) == 32
    blocks = []
    minimum_multipliers = []
    max_coefficient = 0
    for ell in range(1, 32):
        delta = [prefixes[ell][j] - prefixes[ell - 1][j]
                 for j in range(15)]
        assert all(value % 4 == 0 for value in delta)
        coeffs = [sum((delta[j] // 4) * basis[j][k] for j in range(15))
                  for k in range(29)]
        if ell == 31:
            coeffs = [x + y for x, y in zip(coeffs, correction)]
        max_coefficient = max(max_coefficient, max(abs(x) for x in coeffs))
        labels = [sum(coeffs[k] * walks[k]["labels"][r] for k in range(29))
                  for r in range(29)]
        assert labels[0] == (-1 if ell == 31 else 0)
        for r, pair in enumerate(PAIRS, 1):
            if pair in cross:
                expected = 4 * (delta[other.index(pair)] // 4) if pair != ref else 0
            else:
                expected = 4 if ell == 31 and pair in ((0, 1), (2, 3), (4, 5), (6, 7)) else 0
            assert labels[r] == expected
        flow = signed_flow(coeffs)
        needed = 1
        for (_, mask), count in flow.items():
            weight = circ[bin(mask).count("1")]
            assert weight > 0
            while count + needed * weight <= 0:
                needed += 1
        assert needed in (1, 2)
        assert all(count + needed * circ[bin(mask).count("1")] > 0
                   for (_, mask), count in flow.items())
        minimum_multipliers.append(needed)
        blocks.append((coeffs, labels))

    assert len(blocks) == 31
    assert [i + 1 for i, x in enumerate(minimum_multipliers) if x == 2] == [24, 26]
    assert sum(minimum_multipliers) == 33
    assert prefixes[-1] == [0] * 15

    # The signed block corrections telescope; only the old correction
    # remains. Adding K circulations gives the published exact labels.
    aggregate = [sum(block[1][r] for block in blocks) for r in range(29)]
    old = [sum(correction[k] * walks[k]["labels"][r] for k in range(29))
           for r in range(29)]
    assert aggregate == old
    base = data["shortest_terminal_plus_masks"]
    assert len(base) == 11
    state = machine.initial
    for mask in base:
        transitions = {m: nxt for m, nxt, _ in machine.transitions(state)}
        assert mask in transitions
        state = transitions[mask]
    assert state == machine.terminal
    assert all(sum(((mask >> i) & 1) != ((mask >> k) & 1)
                   for mask in base) == 5 for i, k in cross)
    K = 33
    s = 2 + K * data["integer_signed_cycle_correction"]["degree_increment_per_circulation"]
    n = 4 * s + 2
    assert n == len(base) + aggregate[0] + K * data["active_only_total_flow_length"]
    ranks = {row: pos for pos, row in enumerate(machine.kappa_order)}
    for r, (i, k) in enumerate(PAIRS, 1):
        base_distance = sum(((mask >> i) & 1) != ((mask >> k) & 1)
                            for mask in base)
        final_distance = (base_distance + aggregate[r]
                          + K * data["active_only_pair_distance_each"])
        assert final_distance == (2 * s + 1 if (i < 4) != (k < 4)
                                  else 2 * s + 2)
        if (i, k) in cross:
            delta_k = machine.delta[k]
            assert (2 * s + 1 + delta_k) % 2 == 0
            monomial_sign = (-1) ** ((2 * s + 1 + delta_k) // 2)
            prescribed_sign = 1 if ranks[k] > ranks[i] else -1
            assert monomial_sign == prescribed_sign
    print("PASS: 15 integral 4-unit label directions, 31 positive terminal-rooted blocks")
    print(f"PASS: block multipliers sum to {K}; max walk coefficient {max_coefficient}; s={s}, n={n}")


if __name__ == "__main__":
    run()
