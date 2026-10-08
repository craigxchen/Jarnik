"""Cross-check a spread of C++ graph transitions against the Python model.

Usage: python3 this_file.py /path/to/compiled/check_compressed_q2
"""

import json
from collections import Counter
from pathlib import Path
import subprocess
import sys

from check_critical_moment_compressed_automaton import Automaton


def rank_mod_prime(rows, prime=1_000_000_007):
    basis = {}
    for original in rows:
        row = [x % prime for x in original]
        for j in range(len(row)):
            if not row[j]:
                continue
            if j in basis:
                factor = row[j]
                row = [(x-factor*y) % prime for x, y in zip(row, basis[j])]
            else:
                inverse = pow(row[j], -1, prime)
                basis[j] = [(x*inverse) % prime for x in row]
                break
    return len(basis)


def main():
    if len(sys.argv) != 2:
        raise SystemExit("usage: python3 check_critical_moment_compressed_graph_even_q2_reference.py COMPILED_CHECKER")
    fixture = json.loads(Path(__file__).with_name(
        "critical_moment_compressed_graph_even_q2_summary.json"
    ).read_text())
    output = subprocess.run(
        [sys.argv[1], "--sample"], check=True, capture_output=True, text=True
    ).stdout
    machine = Automaton(0, 2)
    lines = iter(output.splitlines())
    state_count = edge_count = 0
    for line in lines:
        values = line.split()
        assert values[0] == "state" and len(values) == 18
        state = (tuple(map(int, values[1:9])), tuple(map(int, values[9:17])))
        count = int(values[17])
        actual = {}
        for _ in range(count):
            edge = next(lines).split()
            assert edge[0] == "edge" and len(edge) == 19
            mask, increments = map(int, edge[1:3])
            destination = (
                tuple(map(int, edge[3:11])),
                tuple(map(int, edge[11:19])),
            )
            assert mask not in actual
            actual[mask] = (destination, increments)
        expected = {
            mask: (destination, sum(bit << i for i, bit in enumerate(distance)))
            for mask, destination, distance in machine.transitions(state)
        }
        assert actual == expected
        assert 0 in actual and 255 in actual
        state_count += 1
        edge_count += count
    assert (state_count, edge_count) == (16, 144)
    full = subprocess.run([sys.argv[1]], check=True, capture_output=True, text=True).stdout
    rows = full.splitlines()
    header = rows[:6]
    assert header[0].split() == [
        "states", str(fixture["reachable_states"]), "edges", str(fixture["labeled_edges"]),
        "processed", str(fixture["reachable_states"]), "closed", "1", "reason", "none",
    ]
    shortest = header[1].split()
    masks = list(map(int, shortest[3:]))
    assert shortest[:3] == ["terminal_depth", str(fixture["shortest_terminal_depth"]), "masks"]
    assert masks == fixture["shortest_terminal_plus_masks"]
    assert header[2].split() == [
        "scc_count", "5", "largest", "322560", "closed_scc", "1", "closed_states", "322560",
    ]
    assert header[3].split() == ["initial_scc_size", "1", "terminal_scc_size", "322560"]
    assert header[4].split() == [
        "all_outdegree9", "1", "all_one_each_size", "1",
        "giant_each_size_permutation", "1", "giant_uniform_label_counts", "1",
    ]
    assert header[5].split() == [
        "active_only_isotropic", "1", "length", str(fixture["active_only_total_flow_length"]),
        "pair_distance", str(fixture["active_only_pair_distance_each"]),
    ]
    state = machine.initial
    distance = [0]*28
    base_flow = Counter()
    for mask in masks:
        transition = {m: (out, inc) for m, out, inc in machine.transitions(state)}
        assert mask in transition
        assert mask not in (0, 255)
        base_flow[(state, mask)] += 1
        state, inc = transition[mask]
        distance = [a+b for a, b in zip(distance, inc)]
    assert state == machine.terminal
    cross = {str(d): sum(d == distance[i] for i, k in enumerate(machine.kind) if k == "C")
             for d in set(distance)}
    assert {k: v for k, v in cross.items() if v} == fixture["shortest_terminal_path_distances"]["cross_pairs"]
    within = {str(d): sum(d == distance[i] for i, k in enumerate(machine.kind) if k != "C")
              for d in set(distance)}
    assert {k: v for k, v in within.items() if v} == fixture["shortest_terminal_path_distances"]["within_pairs"]
    witnesses = []
    for line in rows[6:-1]:
        fields = line.split()
        assert fields[:1] == ["cycle_witness"] and fields[2] == "edge" and fields[4] == "mask" and fields[6] == "labels"
        assert int(fields[1]) == len(witnesses)+1
        path_index = fields.index("path")
        witness = {"edge": int(fields[3]), "mask": int(fields[5]),
                   "labels": list(map(int, fields[7:path_index])),
                   "path_masks": list(map(int, fields[path_index+1:]))}
        assert witness["mask"] not in (0, 255) and len(witness["labels"]) == 29
        assert len(witness["path_masks"]) == witness["labels"][0]
        assert witness["mask"] in witness["path_masks"]
        state = machine.terminal
        label = [0]*29
        for mask in witness["path_masks"]:
            assert mask not in (0, 255)
            transition = {m: (out, inc) for m, out, inc in machine.transitions(state)}
            assert mask in transition
            state, inc = transition[mask]
            label = [label[0]+1]+[x+y for x, y in zip(label[1:], inc)]
        assert state == machine.terminal and label == witness["labels"]
        witnesses.append(witness)
    assert witnesses == fixture["active_cycle_witnesses"]
    assert rows[-1] == "active_cycle_label_rank_mod_1000000007 29"
    assert rank_mod_prime([w["labels"] for w in witnesses]) == 29
    correction = fixture["integer_signed_cycle_correction"]
    coefficients = correction["coefficients"]
    assert len(coefficients) == len(witnesses) == 29
    base_label = [len(masks)]+distance
    corrected_label = [base_label[j]+sum(coefficients[i]*witnesses[i]["labels"][j]
                                      for i in range(29)) for j in range(29)]
    target = [10]+[5 if kind == "C" else 6 for kind in machine.kind]
    assert corrected_label == target == correction["base_path_plus_correction_label"]
    cycle_flow = Counter()
    for coefficient, witness in zip(coefficients, witnesses):
        if coefficient == 0:
            continue
        state = machine.terminal
        for mask in witness["path_masks"]:
            transition = {m: out for m, out, _ in machine.transitions(state)}
            assert mask in transition
            cycle_flow[(state, mask)] += coefficient
            state = transition[mask]
        assert state == machine.terminal
    assert len(cycle_flow) == correction["correction_distinct_directed_edges"]
    assert sum(v < 0 for v in cycle_flow.values()) == correction["correction_negative_directed_edges"]
    assert min(cycle_flow.values()) == correction["correction_minimum_edge_multiplicity"]
    weights = fixture["active_only_integer_flow_per_edge_by_cut_size"]
    combined = dict(cycle_flow)
    for edge, count in base_flow.items():
        combined[edge] = combined.get(edge, 0)+count
    assert correction["minimum_active_circulation_multiplier"] == 1
    assert min(combined[edge]+weights[bin(edge[1]).count("1")] for edge in cycle_flow) == \
        correction["minimum_touched_giant_edge_multiplicity_after_base_and_one_circulation"]
    assert all(combined[edge]+weights[bin(edge[1]).count("1")] > 0 for edge in cycle_flow)
    assert correction["first_result_degree_s"] == \
        correction["base_degree_s"]+correction["degree_increment_per_circulation"]
    assert fixture["active_only_total_flow_length"] == \
        4*correction["degree_increment_per_circulation"]
    assert fixture["active_only_pair_distance_each"] == \
        2*correction["degree_increment_per_circulation"]
    assert correction["first_result_length"] == 4*correction["first_result_degree_s"]+2
    assert correction["first_result_cross_distance"] == 2*correction["first_result_degree_s"]+1
    assert correction["first_result_within_distance"] == 2*correction["first_result_degree_s"]+2
    print("PASS: 16 spread graph states, 144 exact labeled transitions, both constant loops")
    print("PASS: closed graph summary, terminal path and 28 distances, active flow, 29 exact cycle-label witnesses")
    print("PASS: integral s=2 signed correction, K=1 positive active flow, and exact pair labels")


if __name__ == "__main__":
    main()
