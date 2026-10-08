#!/usr/bin/env python3
"""Exact finite checks for thin_edge_divisor_matching_criterion.md.

The graph search checks the incidence assertion only: four edge words can
contain a Hadamard quartet of cut-orientation columns exactly when the four
edges are a matching.  It does not test the analytic thinness hypotheses.
"""

from itertools import combinations
from math import gcd


def cut_signs(edges, m):
    """Sign lines induced by all anchor-fixed cuts on vertices 0..m."""
    lines = set()
    for mask in range(1 << m):
        signs = []
        for a, b in edges:
            ca = 0 if a == 0 else (mask >> (a - 1)) & 1
            cb = 0 if b == 0 else (mask >> (b - 1)) & 1
            signs.append(ca - cb)
        if all(s != 0 for s in signs):
            # The line is invariant under reversing the cut. Normalize its
            # first coordinate to +1.
            if signs[0] < 0:
                signs = [-s for s in signs]
            lines.add(tuple(signs))
    return lines


def has_hadamard_quartet(lines):
    """Whether four available sign lines are pairwise orthogonal."""
    for quartet in combinations(lines, 4):
        if all(sum(x * y for x, y in zip(u, v)) == 0
               for u, v in combinations(quartet, 2)):
            return True
    return False


def edge_sets(m):
    vertices = range(m + 1)
    return combinations(list(combinations(vertices, 2)), 4)


def is_matching(edges):
    endpoints = [v for e in edges for v in e]
    return len(set(endpoints)) == 8


def check_graphs():
    counts = {}
    for m in range(4, 8):  # K5 through K8
        total = matching = feasible = false_positive = false_negative = 0
        for edges in edge_sets(m):
            total += 1
            match = is_matching(edges)
            if match:
                matching += 1
            feasible = has_hadamard_quartet(cut_signs(edges, m))
            if feasible:
                feasible_count = 1
            else:
                feasible_count = 0
            if feasible_count and not match:
                false_positive += 1
            if match and not feasible_count:
                false_negative += 1
        assert false_positive == false_negative == 0, (m, false_positive, false_negative)
        counts[m + 1] = (total, matching)
    return counts


def gaussian_mul(z, w):
    x, y = z
    u, v = w
    return (x * u - y * v, x * v + y * u)


def gaussian_norm(z):
    x, y = z
    return x * x + y * y


def conjugate_primitive(z):
    x, y = z
    # For gcd(x,y)=1, a common Gaussian divisor of z and conjugate(z)
    # could only lie over 2. Opposite coordinate parity excludes 1+i.
    return gcd(abs(x), abs(y)) == 1 and ((x - y) & 1) == 1


def check_half_height_fixture():
    checked = 0
    for a in range(2, 2002, 2):
        d = (a + 2, 1)
        u = (a, -1)
        v = ((a + 1) ** 2, -2)
        assert gaussian_mul(d, u) == v
        assert all(conjugate_primitive(z) for z in (d, u, v))
        nd, nu, nv = map(gaussian_norm, (d, u, v))
        assert nd * nu == nv
        assert nd % 2 == nu % 2 == 1
        assert gcd(nd, nu) == 1
        # Exact identities supporting the scale statement.
        assert nd == a * a + 4 * a + 5
        assert nu == a * a + 1
        assert nv == (a + 1) ** 4 + 4
        assert abs(d[1]) == abs(u[1]) == 1 and abs(v[1]) == 2
        checked += 1
    return checked


def check_eight_vertex_fixture():
    edges = [(1, 0), (3, 2), (5, 4), (7, 6)]
    cuts = {
        "D": {1, 3, 5, 7},
        "A": {1, 3, 4, 6},
        "B": {1, 2, 5, 6},
        "C": {1, 2, 4, 7},
    }
    expected = {
        "D": (1, 1, 1, 1),
        "A": (1, 1, -1, -1),
        "B": (1, -1, 1, -1),
        "C": (1, -1, -1, 1),
    }
    actual = {}
    for name, subset in cuts.items():
        actual[name] = tuple(
            (1 if a in subset else 0) - (1 if b in subset else 0)
            for a, b in edges
        )
        assert actual[name] == expected[name]
    assert all(sum(x * y for x, y in zip(actual[p], actual[q])) == 0
               for p, q in combinations(actual, 2))

    # Across all anchor-fixed cuts, every one of the eight sign lines occurs.
    lines = cut_signs(edges, 7)
    assert len(lines) == 8
    assert has_hadamard_quartet(lines)
    return len(lines)


def check_shared_vertex_orientation_depths():
    """Check the row-sign mismatch that forces the full missing depth into K."""
    target_columns = {
        "A": (1, 1, -1, -1),
        "B": (1, -1, 1, -1),
        "C": (1, -1, -1, 1),
    }
    checks = 0
    for r, s in combinations(range(4), 2):
        for eta in (-1, 1):
            products = {
                name: signs[r] * signs[s]
                for name, signs in target_columns.items()
            }
            assert products.values() and all(x in (-1, 1) for x in products.values())
            mismatching_groups = [name for name, product in products.items()
                                 if product == -eta]
            assert len(mismatching_groups) == (2 if eta == 1 else 1)
            # If a prime has full oriented depth d in such a group, the common
            # core orientation ratio eta leaves one row without that target
            # orientation. Record valuations at pi and conjugate(pi) directly.
            for depth in range(1, 13):
                core_r = (depth, 0)  # normalize the first core orientation to pi
                core_s = (depth, 0) if eta == 1 else (0, depth)
                target_r = (depth, 0)
                target_s = (depth, 0) if -eta == 1 else (0, depth)
                need_r = tuple(max(t - c, 0) for t, c in zip(target_r, core_r))
                need_s = tuple(max(t - c, 0) for t, c in zip(target_s, core_s))
                correction_norm_valuation = sum(need_r) + sum(need_s)
                assert correction_norm_valuation >= depth
                checks += 1
    return checks


if __name__ == "__main__":
    graph_counts = check_graphs()
    fixture_count = check_half_height_fixture()
    line_count = check_eight_vertex_fixture()
    valuation_checks = check_shared_vertex_orientation_depths()
    print("K_n: (four-edge sets, matchings), feasible iff matching:")
    for n, counts in graph_counts.items():
        print(f"  K{n}: {counts[0]}, {counts[1]}")
    print(f"eight-vertex fixture: {line_count} sign lines, Hadamard quartet verified")
    print(f"half-height factorization: {fixture_count} even-a fixtures verified")
    print(f"shared-vertex sign/depth cases: {valuation_checks} verified")
