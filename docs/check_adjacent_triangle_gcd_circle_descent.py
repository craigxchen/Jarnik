#!/usr/bin/env python3
"""Exact checks for adjacent_triangle_gcd_circle_descent_obstruction.md."""

from __future__ import annotations

from functools import reduce
from itertools import combinations
from math import gcd

from check_complement_reflection_invariant_relations import (
    exact_div, gaussian_gcd, norm,
)
from check_five_point_affine_shape_cubic_family import (
    SUBSETS, evaluate, forward_point, reverse_degree_six, triangle,
)


def extended_gcd_list(values: list[int]) -> tuple[int, list[int]]:
    g = 0
    coefficients: list[int] = []
    for value in values:
        r0, r1 = g, value
        old_s, s = 1, 0
        old_t, t = 0, 1
        while r1:
            q = r0 // r1
            r0, r1 = r1, r0-q*r1
            old_s, s = s, old_s-q*s
            old_t, t = t, old_t-q*t
        if r0 < 0:
            r0, old_s, old_t = -r0, -old_s, -old_t
        coefficients = [old_s*c for c in coefficients]+[old_t]
        g = r0
    assert sum(c*v for c, v in zip(coefficients, values)) == g
    return g, coefficients


def lattice_basis(vectors: list[tuple[int, int]], index: int) -> tuple[tuple[int, int], tuple[int, int]]:
    x_gcd, coefficients = extended_gcd_list([v[0] for v in vectors])
    assert x_gcd > 0 and index % x_gcd == 0
    y0 = sum(c*v[1] for c, v in zip(coefficients, vectors))
    vertical = index//x_gcd
    y0 %= vertical
    first, second = (x_gcd, y0), (0, vertical)
    for x, y in vectors:
        assert x % x_gcd == 0
        assert (x_gcd*y-y0*x) % index == 0
    return first, second


def main() -> None:
    reverse = [tuple(map(reverse_degree_six, forward_point(s))) for s in SUBSETS]
    checked = 0
    for k in list(range(1, 41))+[100, 1000]:
        T = 600*k
        raw = [(evaluate(x, T), evaluate(y, T)) for x, y in reverse]
        points = [exact_div(z, (720, 0)) for z in raw]
        assert norm(reduce(gaussian_gcd, points)) == 1
        assert len({norm(z) for z in points}) == 1

        triples = list(combinations(range(5), 3))
        determinants = [abs(triangle(points, triple)) for triple in triples]
        full_index = reduce(gcd, determinants)
        assert full_index == T//15
        adjacent = [abs(triangle(points, triple))
                    for triple in ((0, 1, 2), (1, 2, 3), (2, 3, 4))]
        assert reduce(gcd, adjacent) == full_index
        primitive = [value//full_index for value in adjacent]
        expected = [3*(T*T+4)//4, 7*T*T//12, 25*(T*T+36)//36]
        assert primitive == expected
        assert reduce(gcd, primitive) == 1
        assert 3*max(primitive) < 4*min(primitive)

        differences = [(z[0]-points[0][0], z[1]-points[0][1]) for z in points[1:]]
        all_minors = [abs(a[0]*b[1]-a[1]*b[0]) for a, b in combinations(differences, 2)]
        assert reduce(gcd, all_minors) == full_index
        col1, col2 = lattice_basis(differences, full_index)
        metric = (
            (norm(col1), col1[0]*col2[0]+col1[1]*col2[1]),
            (col1[0]*col2[0]+col1[1]*col2[1], norm(col2)),
        )
        assert metric[0][0]*metric[1][1]-metric[0][1]**2 == full_index**2
        # If this basis were a similarity and its center integral, it would
        # exhibit a common Gaussian divisor, contrary to the checked gcd.
        similarity = metric[0][1] == 0 and metric[0][0] == metric[1][1]
        det_a = col1[0]*col2[1]-col1[1]*col2[0]
        center_integral = (
            (col2[1]*points[0][0]-col2[0]*points[0][1]) % det_a == 0
            and (-col1[1]*points[0][0]+col1[0]*points[0][1]) % det_a == 0
        )
        assert not (similarity and center_integral)

        x = T*T
        aa = 4*(x+25)*(x+36)
        bb = 10*(x*x+43*x+360)
        cc = 25*(x+9)*(x+16)
        assert aa*cc-bb*bb == (180*T)**2
        conic_content = gcd(gcd(aa, 2*bb), cc)
        assert conic_content == 3600
        assert max(aa, 2*bb, cc)//conic_content == cc//3600
        basis = [(points[j][0]-points[0][0], points[j][1]-points[0][1])
                 for j in (1, 2)]
        assert norm(basis[0])*3600 == (x+4)*aa
        assert norm(basis[1])*3600 == (x+4)*cc
        assert (basis[0][0]*basis[1][0]+basis[0][1]*basis[1][1])*3600 == (x+4)*bb
        checked += 1

    print(f"PASS: {checked} primitive endpoint family members with adjacent gcd T/15 and ratio <4/3.")
    print("PASS: full difference-lattice index, determinant-g^2 metric, and conic scale identities.")


if __name__ == "__main__":
    main()
