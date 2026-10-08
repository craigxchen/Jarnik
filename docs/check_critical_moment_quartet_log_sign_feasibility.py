"""Exact 29-cycle Monge label and positive-flow audit; no huge Euler word."""

from collections import Counter
from itertools import combinations
import json
from pathlib import Path

from check_critical_moment_compressed_automaton import Automaton, PAIRS


MONGE_CYCLE_COEFFICIENTS = (
    0, -2, 25, -19, 0, -18, 14, -20, 23, 0, 10, 6, -9, -2, 1,
    -3, 3, 3, -2, -3, 1, 0, -3, -3, -3, -3, 1, 3, 0,
)


def run():
    data = json.loads((Path(__file__).with_name(
        "critical_moment_compressed_graph_even_q2_summary.json")).read_text())
    machine = Automaton(0, 2)
    cycles = data["active_cycle_witnesses"]
    assert len(cycles) == len(MONGE_CYCLE_COEFFICIENTS) == 29
    correction = data["integer_signed_cycle_correction"]["coefficients"]
    assert len(correction) == 29
    circular = data["active_only_integer_flow_per_edge_by_cut_size"]
    assert circular == [0, 10, 28, 56, 70, 56, 28, 10, 0]

    target = [0] + [(-4 * i * (k - 4) if i < 4 <= k else 0)
                    for i, k in PAIRS]
    labels = [sum(c * cycle["labels"][j]
                  for c, cycle in zip(MONGE_CYCLE_COEFFICIENTS, cycles))
              for j in range(29)]
    assert labels == target

    cycle_edges = []
    for cycle in cycles:
        state = machine.terminal
        edges = []
        replayed_labels = [0] * 29
        for mask in cycle["path_masks"]:
            next_states = {m: (nxt, inc) for m, nxt, inc in machine.transitions(state)}
            assert mask in next_states and mask not in (0, 255)
            edges.append((state, mask))
            state, increments = next_states[mask]
            replayed_labels[0] += 1
            for j, inc in enumerate(increments, 1):
                replayed_labels[j] += inc
        assert state == machine.terminal
        assert replayed_labels == cycle["labels"]
        cycle_edges.append(edges)

    def edge_correction(coefficients):
        result = Counter()
        for c, edges in zip(coefficients, cycle_edges):
            for edge in edges:
                result[edge] += c
        return result

    first = edge_correction(MONGE_CYCLE_COEFFICIENTS)
    second = edge_correction(tuple(a - b for a, b in
                                   zip(correction, MONGE_CYCLE_COEFFICIENTS)))
    for multiplier, edge_weights, expected_min in ((4, first, 2), (3, second, 3)):
        touched = [weight + multiplier * circular[bin(mask).count("1")]
                   for (_, mask), weight in edge_weights.items() if weight]
        assert min(touched) == expected_min
        assert all(weight > 0 for weight in touched)

    ranks = {row: pos for pos, row in enumerate(machine.kappa_order)}
    base_masks = data["shortest_terminal_plus_masks"]
    assert len(base_masks) == 11
    base_cross = [sum(((mask >> i) & 1) != ((mask >> k) & 1)
                      for mask in base_masks) for i in range(4) for k in range(4, 8)]
    assert base_cross == [5] * 16
    for i, j in combinations(range(4), 2):
        for k, l in combinations(range(4, 8), 2):
            c_i, c_j, c_k, c_l = (ranks[r] for r in (i, j, k, l))
            w = (c_i - c_j) * (c_k - c_l)
            v = (c_i - c_l) * (c_j - c_k)
            assert w * v > 0
            d = {(r, b): target[1 + PAIRS.index((r, b))] for r in (i, j)
                 for b in (k, l)}
            h = (d[i, l] + d[j, k] - d[i, k] - d[j, l]) // 2
            assert h == 2 * (j - i) * (l - k) > 0

    s = 2 + 7 * data["integer_signed_cycle_correction"]["degree_increment_per_circulation"]
    n = 4 * s + 2
    assert n == len(base_masks) + 7 * data["active_only_total_flow_length"] - 1
    print(f"PASS: integer Monge labels, two positive giant flows, 36 positive quartet prefix counts; s={s}, n={n}.")


if __name__ == "__main__":
    run()
