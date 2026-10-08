"""Boundary-node restriction ranks for the five-row coefficient lemma.

The arithmetic modulus is proved in the companion note. This checks the
full invariant space and both incident restrictions at all 15 nodes.
"""

from itertools import combinations

from check_odd_support_relation_rigidity import graphs, polynomial_rank, rank
from check_smaller_cut_relation_rank import restrict_graph


def linear_coefficients(poly):
    assert all(sum(monomial) == 1 for monomial in poly)
    result = [poly.get(tuple(int(i == j) for i in range(5)), 0)
              for j in range(5)]
    assert sum(result) == 0
    return result


def main():
    monomials = list(graphs((2,) * 5))
    assert len(monomials) == 22
    assert polynomial_rank([restrict_graph(edges, set(), 5)
                            for edges in monomials]) == 6
    boundaries = list(combinations(range(5), 2))
    restrictions = {
        cut: [linear_coefficients(restrict_graph(edges, set(cut), 5))
              for edges in monomials]
        for cut in boundaries
    }
    nodes = 0
    for first, second in combinations(boundaries, 2):
        if not set(first).isdisjoint(second):
            continue
        remaining = next(i for i in range(5) if i not in first + second)
        left, right = restrictions[first], restrictions[second]
        # At the node, outside second-pair coordinates are zero and the
        # final coordinate is one. The opposite boundary chart gives
        # the same node functional up to the SL2 swap sign.
        node_values = [row[remaining] for row in left]
        assert any(node_values)
        assert [row[remaining] for row in right] == [-x for x in node_values]
        anchor = next(j for j, value in enumerate(node_values) if value)
        left_slopes, right_slopes, paired = [], [], []
        for j in range(len(monomials)):
            if j == anchor:
                continue
            # These vectors span the kernel of the node-value functional.
            l = [node_values[anchor] * a - node_values[j] * b
                 for a, b in zip(left[j], left[anchor])]
            r = [node_values[anchor] * a - node_values[j] * b
                 for a, b in zip(right[j], right[anchor])]
            assert l[remaining] == r[remaining] == 0
            assert all(l[i] == 0 for i in first)
            assert all(r[i] == 0 for i in second)
            assert l[second[0]] == -l[second[1]]
            assert r[first[0]] == -r[first[1]]
            left_slopes.append([l[second[0]]])
            right_slopes.append([r[first[0]]])
            paired.append([l[second[0]], r[first[0]]])
        assert rank(left_slopes) == rank(right_slopes) == 1
        assert rank(paired) == 2
        nodes += 1
    assert nodes == 15
    assert 6 + 2 + 2 == 10
    assert 10 + (10 + 10 * 9) == 110
    print("Passed: all 22 graph coordinates span the six-dimensional invariant space.")
    print("Passed: all 15 nodes give single-difference restrictions on both incident lines.")
    print("Passed: each individual slope has rank one; both slopes together have rank two.")


if __name__ == "__main__":
    main()
