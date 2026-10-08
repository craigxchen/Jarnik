#!/usr/bin/env python3
"""Exact constants for the Chebyshev relaxation-feasibility theorem.

Checks symbolic numerator coefficients and saved mask multiplicities.
Does not construct or assert any solution of the odd-moment equations.
"""

from collections import Counter
from fractions import Fraction as F
from itertools import combinations
from pathlib import Path
import json

from check_critical_moment_multirow_christoffel import (
    determinant, inverse, multiply, transpose,
)


def check_integer_gram_identity():
    """Audit the determinant formula, including zero-support obstructions."""
    s, n = 3, 14
    walsh = [[(-1) ** bin(i & j).count('1') for j in range(16)]
             for i in range(8)]
    v = [row[2:] for row in walsh]
    b = [[v[i][j]-v[0][j] for j in range(n)] for i in range(1, 8)]
    gram = multiply(b, transpose(b))
    det_gram = determinant(gram)
    assert det_gram == 8*(4*s+4)**6*(4*s-4)
    q = multiply(multiply(transpose(b), inverse(gram)), b)
    zeros = positives = 0
    for chosen in combinations(range(n), s):
        complement = [j for j in range(n) if j not in chosen]
        restricted = [[row[j] for j in complement] for row in b]
        integer_numerator = determinant(multiply(restricted, transpose(restricted)))
        principal = [[F(j == k)-q[j][k] for k in chosen] for j in chosen]
        value = determinant(principal)
        assert value == integer_numerator/det_gram
        assert integer_numerator.denominator == 1 and integer_numerator >= 0
        if integer_numerator:
            positives += 1
            assert value >= 1/det_gram
        else:
            zeros += 1
    assert zeros > 0 and positives > 0 and zeros+positives == 364
    return zeros, positives


def main():
    zero_minors, positive_minors = check_integer_gram_identity()
    # Scalar margin at s=5+t: numerator 34t^2+213t+80.
    coefficients = [34, 2 * 34 * 5 - 127, 34 * 25 - 127 * 5 - 135]
    assert coefficients == [34, 213, 80]
    assert all(c > 0 for c in coefficients)
    assert F(21, 38) < F(3, 5)

    # Exponential versus determinant denominator, then monotone propagation.
    s = 100
    determinant = 8 * (4 * s + 4) ** 6 * (4 * s - 4)
    assert 3 ** s * determinant < 5 ** s
    assert 3 * 100 ** 7 < 5 * 99 ** 7

    path = Path(__file__).with_name('critical_moment_compressed_graph_even_q2_summary.json')
    data = json.loads(path.read_text())
    correction = Counter(data['shortest_terminal_plus_masks'])
    coefficients = data['integer_signed_cycle_correction']['coefficients']
    witnesses = data['active_cycle_witnesses']
    assert len(coefficients) == len(witnesses) == 29
    for coefficient, witness in zip(coefficients, witnesses):
        for mask in witness['path_masks']:
            correction[mask] += coefficient
    assert sum(correction.values()) == 10
    assert min(correction.values()) == -6
    assert max(correction.values()) == 5
    assert correction[0] == correction[255] == 0

    circulation_per_mask = []
    for count, weight in zip(data['edges_per_specific_label_in_closed_scc'],
                             data['active_only_integer_flow_per_edge_by_cut_size']):
        circulation_per_mask.append(count * weight)
    assert circulation_per_mask == [0, 403200, 322560, 322560, 322560,
                                    322560, 322560, 403200, 0]
    assert data['integer_signed_cycle_correction']['degree_increment_per_circulation'] == 20805120

    # rho<39/40 iff 39^2(4s-4)^2 - 40^2*768*s*m > 0.
    # Let k=1+t, so s=A*t+A+2 and m=B*t+B+5.
    a, b = 20805120, 403200
    c, d = 4 * a, 4 * a + 4
    quadratic = [
        39 ** 2 * c ** 2 - 40 ** 2 * 768 * a * b,
        39 ** 2 * 2 * c * d
        - 40 ** 2 * 768 * (a * (b + 5) + b * (a + 2)),
        39 ** 2 * d ** 2 - 40 ** 2 * 768 * (a + 2) * (b + 5),
    ]
    assert all(value > 0 for value in quadratic)
    assert 3 ** 64 * 40 ** 7 < 5 ** 64

    print('PASS: all scalar bounds for s>=5; all determinant bounds under')
    print('      the support condition for s>=100; exact mask correction;')
    print('      rho<39/40 for every established degree s=2+20805120k, k>=1.')
    print('rho numerator coefficients at k=1+t:', quadratic)
    print('Integer Gram identity: {} zero and {} positive subset minors.'.format(
        zero_minors, positive_minors))
    print('SCOPE: relaxations only; no original moment-equation feasibility claim.')


if __name__ == '__main__':
    main()
