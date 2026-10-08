"""Exact exterior-coordinate diagnostic for the SL_2, d=3 kernel.

Enumerates polynomial coordinates, with integer arithmetic and no residue
specialization. Also verifies the stated fusion-budget normalizations.
"""

from collections import Counter
from itertools import combinations, permutations

from check_higher_rank_terminal_fusion_grid import case


def popcount(mask):
    return bin(mask).count("1")


def noncrossing_matchings(labels):
    if not labels:
        yield ()
        return
    for j in range(1, len(labels), 2):
        for inside in noncrossing_matchings(labels[1:j]):
            for outside in noncrossing_matchings(labels[j + 1:]):
                yield ((labels[0], labels[j]),) + inside + outside


def matching_polynomial(matching):
    """Coefficients of product(z_j-z_i), indexed by squarefree masks."""
    polynomial = {0: 1}
    for i, j in matching:
        polynomial = {
            mask | (1 << k): coefficient * sign
            for mask, coefficient in polynomial.items()
            for k, sign in ((j, 1), (i, -1))
        }
    return polynomial


SIGNED_PERMUTATIONS = tuple(
    (p, (-1) ** sum(p[i] > p[j] for i in range(4)
                   for j in range(i + 1, 4)))
    for p in permutations(range(4))
)


def determinant_four(matrix):
    result = 0
    for permutation, sign in SIGNED_PERMUTATIONS:
        term = sign
        for i in range(4):
            term *= matrix[i][permutation[i]]
        result += term
    return result


def main():
    matchings = tuple(noncrossing_matchings(tuple(range(6))))
    assert len(matchings) == 5
    polynomials = tuple(matching_polynomial(item) for item in matchings)
    masks = tuple(sum(1 << i for i in subset)
                  for subset in combinations(range(6), 3))
    embedding = tuple(tuple(poly.get(mask, 0) for poly in polynomials)
                      for mask in masks)
    coalitions = tuple(range(1 << 6))
    common_orders = tuple(max(0, popcount(t) - 3) for t in coalitions)
    assert sum(common_orders) == 30

    # Verify these are the actual common polynomial orders of the tuple.
    for t, expected in zip(coalitions, common_orders):
        assert min(popcount(mask & t)
                   for poly in polynomials for mask in poly) == expected

    distribution = Counter()
    zero_count = 0
    witness_masks = (7, 11, 21, 25)
    witness_seen = False
    for selected in combinations(range(20), 4):
        # Exterior projection of the signed dual of the evaluation row.
        multipliers = tuple(
            (-1) ** omitted * determinant_four([
                [embedding[row][column] for column in range(5)
                 if column != omitted] for row in selected
            ]) for omitted in range(5)
        )
        coefficients = tuple(
            sum(multipliers[j] * embedding[i][j] for j in range(5))
            for i in range(20)
        )
        support = tuple(mask for mask, coefficient in zip(masks, coefficients)
                        if coefficient)
        if not support:
            zero_count += 1
            continue
        orders = tuple(min(popcount(mask & t) for mask in support)
                       for t in coalitions)
        residual = tuple(order - common for order, common
                         in zip(orders, common_orders))
        assert min(residual) >= 0
        assert sum(orders) <= 48
        distribution[sum(residual)] += 1
        if tuple(masks[i] for i in selected) == witness_masks:
            witness_seen = True
            assert multipliers == (0, -1, 1, 0, 0)
            positive = set(witness_masks) | {63 ^ mask for mask in witness_masks}
            assert residual == tuple(int(t in positive) for t in coalitions)

    assert witness_seen
    assert zero_count == 1725
    assert distribution == Counter({18: 2880, 8: 240})
    assert zero_count + sum(distribution.values()) == 4845
    print("(2,3): 4845 raw projections; 1725 zero; "
          "240 have sum r=8, 2880 have sum r=18.")

    from math import comb
    expected = {
        (2, 3): (4, 1, 48, 30, 18),
        (3, 3): (21, 8, 4608, 2520, 2088),
        (3, 4): (341, 55, 337920, 166320, 171600),
    }
    for parameters, values in expected.items():
        result = case(*parameters)
        m, delta = result["m"], result["delta"]
        raw = delta * m * (1 << (m - 3))
        content = sum(comb(m, s) * value
                      for s, value in enumerate(result["nu"]))
        actual = result["h"], delta, raw, content, result["B"]
        assert actual == values
        print("{}: h={}, delta={}, raw={}, fusion content={}, B={}".format(
            parameters, *actual))


if __name__ == "__main__":
    main()
