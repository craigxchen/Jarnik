"""Exact Petersen cycle identities and tree propagation certificates."""

from fractions import Fraction

from check_five_row_boundary_slope_relations import EXPECTED


VERTICES = ("12", "13", "14", "15", "23", "24", "25", "34", "35", "45")
TREE = (("12", "34"), ("12", "35"), ("12", "45"),
        ("13", "24"), ("13", "25"), ("13", "45"),
        ("14", "23"), ("14", "25"), ("15", "23"))
CHORDS = (("14", "35"), ("15", "24"), ("15", "34"),
          ("23", "45"), ("24", "35"), ("25", "34"))


def add(left, right):
    result = dict(left)
    for monomial, coefficient in right.items():
        result[monomial] = result.get(monomial, 0) + coefficient
    return {m: c for m, c in result.items() if c}


def scale(poly, coefficient):
    return {m: coefficient * c for m, c in poly.items() if coefficient * c}


def mul(left, right):
    result = {}
    for a, ca in left.items():
        for b, cb in right.items():
            monomial = tuple(x + y for x, y in zip(a, b))
            result[monomial] = result.get(monomial, 0) + ca * cb
    return {m: c for m, c in result.items() if c}


def form(vector):
    return {
        tuple(int(i == j) for i in range(6)): coefficient
        for j, coefficient in enumerate(vector) if coefficient
    }


def slope_forms(label):
    pair = tuple(int(character) - 1 for character in label)
    alpha, beta = EXPECTED[pair]
    return form(alpha), form(beta), add(form(alpha), form(beta))


def node_position(vertex, neighbor):
    assert set(vertex).isdisjoint(neighbor)
    remaining = next(iter(set("12345") - set(vertex) - set(neighbor)))
    return sorted(set("12345") - set(vertex)).index(remaining)


def node_polynomial(vertex, neighbor):
    alpha, beta, total = slope_forms(vertex)
    return (alpha, beta, scale(total, -1))[node_position(vertex, neighbor)]


def tree_path(first, last):
    queue = [(first,)]
    seen = {first}
    while queue:
        path = queue.pop(0)
        if path[-1] == last:
            return path
        for left, right in TREE:
            neighbor = right if left == path[-1] else left if right == path[-1] else None
            if neighbor is not None and neighbor not in seen:
                seen.add(neighbor)
                queue.append(path + (neighbor,))
    raise AssertionError("Tree is disconnected")


def check_tree_and_cycles():
    graph_edges = {frozenset((a, b)) for i, a in enumerate(VERTICES)
                   for b in VERTICES[i+1:] if set(a).isdisjoint(b)}
    assert {frozenset(edge) for edge in TREE + CHORDS} == graph_edges
    assert len(TREE) == 9 and len(graph_edges) == 15
    assert all(tree_path("12", vertex) for vertex in VERTICES)
    lengths = []
    for first, last in CHORDS:
        cycle = tree_path(first, last)
        lengths.append(len(cycle))
        outgoing = incoming = {(0,) * 6: 1}
        for index, vertex in enumerate(cycle):
            outgoing = mul(outgoing, node_polynomial(vertex, cycle[(index+1) % len(cycle)]))
            incoming = mul(incoming, node_polynomial(vertex, cycle[index-1]))
        assert add(outgoing, scale(incoming, -(-1)**len(cycle))) == {}
    assert lengths == [6, 6, 8, 5, 5, 5]

    # Tree propagation on arbitrary nonnodal local slope triples. Even
    # when individual chord gluings fail, their defects sum to zero.
    # Thus any five zero defects force the sixth to be zero.
    for offset in (2, 7, 13):
        local = {vertex: (Fraction(offset+j), Fraction(1), Fraction(-offset-j-1))
                 for j, vertex in enumerate(VERTICES)}
        scales = {"12": Fraction(1)}
        while len(scales) < len(VERTICES):
            for a, b in TREE:
                if (a in scales) == (b in scales):
                    continue
                if b in scales:
                    a, b = b, a
                scales[b] = -scales[a]*local[a][node_position(a, b)]/local[b][node_position(b, a)]
        def defect(a, b):
            return scales[a]*local[a][node_position(a, b)] + scales[b]*local[b][node_position(b, a)]
        assert all(defect(a, b) == 0 for a, b in TREE)
        assert sum(defect(a, b) for a, b in CHORDS) == 0


def main():
    check_tree_and_cycles()
    alpha13, beta13, _ = slope_forms("13")
    alpha14, beta14, sum14 = slope_forms("14")
    alpha23, beta23, sum23 = slope_forms("23")
    alpha25, beta25, sum25 = slope_forms("25")
    alpha45, beta45, _ = slope_forms("45")

    # Clear the denominators in
    # r13*r45*(r23+1)*(r25+1)-r23*(r14+1).
    left = mul(mul(mul(alpha13, alpha45), sum23), mul(sum25, beta14))
    right = mul(mul(mul(alpha23, sum14), beta13), mul(beta45, beta25))
    assert add(left, scale(right, -1)) == {}

    # The five-cycle is nonzero for this concrete section, so the ratio form
    # is legitimate as well as its cleared polynomial form.
    q = (1, 2, 3, 4, 5, 6)
    values = {}
    for label in ("13", "14", "23", "25", "45"):
        alpha, beta, total = slope_forms(label)
        # Evaluate directly using the linear-form support.
        av = sum(coefficient * q[index]
                 for monomial, coefficient in alpha.items()
                 for index, exponent in enumerate(monomial) if exponent)
        bv = sum(coefficient * q[index]
                 for monomial, coefficient in beta.items()
                 for index, exponent in enumerate(monomial) if exponent)
        tv = av + bv
        assert av and bv and tv
        values[label] = (av, bv)
    r = {label: Fraction(alpha, beta) for label, (alpha, beta) in values.items()}
    assert r["13"] * r["45"] * (r["23"] + 1) * (r["25"] + 1) \
        - r["23"] * (r["14"] + 1) == 0
    print("PASS: the 23-14-25-13-45 five-cycle identity holds symbolically and numerically.")
    print("PASS: spanning tree, six exact fundamental cycle identities and their stated lengths.")
    print("PASS: exact rational tree propagation and conservation of the six chord defects.")


if __name__ == "__main__":
    main()
