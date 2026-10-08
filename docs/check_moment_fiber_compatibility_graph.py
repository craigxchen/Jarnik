"""Exhaustive finite compatibility-graph bounds for odd-moment fibers.

The graph is a necessary-condition relaxation: actual fibers must be
cliques, while supports larger than 2*s+2 are intentionally unrestricted.
"""

from itertools import product


def pattern_critical(s):
    """The length 2*s+1 pattern ++--++--..., up to negation."""
    return tuple(1 if (i // 2) % 2 == 0 else -1 for i in range(2 * s + 1))


def pattern_successor_a(s):
    """The length 2*s+2 pattern +--++--..., up to negation."""
    t = 2 * s + 2
    return tuple(
        1 if i == 0 else (-1 if ((i - 1) // 2) % 2 == 0 else 1)
        for i in range(t)
    )


def pattern_successor_b(s):
    """The length 2*s+2 pattern +++--++--..., up to negation."""
    t = 2 * s + 2
    return tuple(
        1 if i < 3 else (-1 if ((i - 3) // 2) % 2 == 0 else 1)
        for i in range(t)
    )


def compatible(a, b, s):
    support = [i for i, (x, y) in enumerate(zip(a, b)) if x != y]
    t = len(support)
    if t <= 2 * s:
        return False
    if t > 2 * s + 2:
        return True

    observed = tuple(a[i] for i in support)
    if t == 2 * s + 1:
        pattern = pattern_critical(s)
        return observed in (pattern, tuple(-x for x in pattern))
    if t == 2 * s + 2:
        patterns = (pattern_successor_a(s), pattern_successor_b(s))
        return any(observed in (p, tuple(-x for x in p)) for p in patterns)
    raise AssertionError("unreachable support size")


def build_graph(m, s):
    vertices = list(product((-1, 1), repeat=m))
    n = len(vertices)
    adjacency = [0] * n
    for i, a in enumerate(vertices):
        mask = 0
        for j, b in enumerate(vertices):
            if i != j and compatible(a, b, s):
                mask |= 1 << j
        adjacency[i] = mask
    return vertices, adjacency


def color_order(candidates, adjacency):
    """Greedy coloring upper bound, returned in search order."""
    remaining = candidates
    order = []
    bounds = []
    color = 0
    while remaining:
        color += 1
        available = remaining
        while available:
            bit = available & -available
            vertex = bit.bit_length() - 1
            available ^= bit
            remaining ^= bit
            order.append(vertex)
            bounds.append(color)
            available &= ~adjacency[vertex]
    return order, bounds


def maximum_clique(adjacency):
    """Exact maximum clique by coloring branch and bound."""
    best = []
    nodes = 0

    def expand(clique, candidates):
        nonlocal best, nodes
        nodes += 1
        if len(clique) > len(best):
            best = clique[:]
        if not candidates:
            return

        order, bounds = color_order(candidates, adjacency)
        for vertex, bound in reversed(list(zip(order, bounds))):
            if len(clique) + bound <= len(best):
                return
            bit = 1 << vertex
            if not (candidates & bit):
                continue
            expand(clique + [vertex], candidates & adjacency[vertex])
            candidates ^= bit

    expand([], (1 << len(adjacency)) - 1)
    return best, nodes


def check_case(m, s, expected):
    vertices, adjacency = build_graph(m, s)
    clique, nodes = maximum_clique(adjacency)
    assert len(clique) == expected, (m, s, len(clique), expected)
    assert len(set(clique)) == len(clique)
    for i, u in enumerate(clique):
        for v in clique[i + 1 :]:
            assert (adjacency[u] >> v) & 1
    print(
        f"PASS: G({m},{s}) has maximum clique {expected}; "
        f"{len(vertices)} vertices, {nodes} search nodes."
    )


check_case(6, 1, 5)
check_case(10, 2, 7)
check_case(8, 2, 3)
