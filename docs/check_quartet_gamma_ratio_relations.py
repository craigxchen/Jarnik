#!/usr/bin/env python3
"""Small exact fixtures for ratios of quartet Gaussian contents.

The gcd gamma_Q is only defined up to a unit, independently for each quartet.
Consequently raw additive identities and complex phases are not invariants.
This checker displays one deterministic Euclidean-GCD choice and tests the
unit-invariant weaker condition: two ratios can be made parallel by units iff
one component of conjugate(a)*b vanishes.
"""

from __future__ import annotations

from itertools import combinations
from math import atan2, gcd

from check_four_vertical_lcm_search import exact_div, gcd_gaussian, mul, norm
from check_moving_quadrilateral_gamma_content import (primitive_tuple,
                                                       quartet_gamma)

Gaussian = tuple[int, int]


def conj(z: Gaussian) -> Gaussian:
    return (z[0], -z[1])


def sub(a: Gaussian, b: Gaussian) -> Gaussian:
    return (a[0] - b[0], a[1] - b[1])


def canonical(z: Gaussian) -> Gaussian:
    """Choose one associate, solely to make printed fixtures reproducible."""
    associates = (z, (-z[1], z[0]), (-z[0], -z[1]), (z[1], -z[0]))
    return min(associates, key=lambda q: (q[0] < 0, q[1] < 0, q[0], q[1]))


def product(a: Gaussian, b: Gaussian) -> Gaussian:
    return mul(a, b)


def fixture(name: str, raw: tuple[Gaussian, ...]) -> None:
    points = tuple(sorted(primitive_tuple(raw), key=lambda q: atan2(q[1], q[0])))
    gammas = []
    for inds in combinations(range(len(points)), 4):
        value, gamma, _, _ = quartet_gamma(tuple(points[i] for i in inds))
        gammas.append((inds, value, gamma))
    G = gammas[0][2]
    for _, _, gamma in gammas[1:]:
        G = gcd_gaussian(G, gamma)
    ratios = [(inds, value, exact_div(gamma, G))
              for inds, value, gamma in gammas]
    ratios = [(inds, value, canonical(ratio)) for inds, value, ratio in ratios]
    determinants = []
    unit_parallel = 0
    for a, b in combinations(ratios, 2):
        c = product(conj(a[2]), b[2])
        determinants.append(abs(c[1]))
        unit_parallel += int(c[0] == 0 or c[1] == 0)
    nonzero = [d for d in determinants if d]
    print(f"{name}: points={points}, G={G}, NormG={norm(G)}")
    print(f"  ratio norms={sorted(value // norm(G) for _, value, _ in ratios)}")
    print(f"  normalized ratios={[r for _, _, r in ratios]}")
    print(f"  |Im(conj(a)b)| range=({min(determinants)},{max(determinants)}), "
          f"min_nonzero={min(nonzero) if nonzero else 0}, "
          f"unit-parallel-pairs={unit_parallel}/{len(determinants)}")


def pi_power(e: int) -> tuple[Gaussian, ...]:
    z = (1, 0)
    for _ in range(e):
        z = mul(z, (2, 1))
    return (z, (-z[1], z[0]), (-z[0], -z[1]), (z[1], -z[0]),
            (z[0], -z[1]))


def main() -> None:
    fixture("pi^1 five", pi_power(1))
    fixture("pi^2 five", pi_power(2))
    fixture("pi^3 five", pi_power(3))
    # Six-point fixture from the same complete circle, retaining one global
    # normalization for all six rows.
    fixture("N=125 six", ((-11, -2), (-10, -5), (-5, -10), (-2, -11),
                           (2, -11), (5, -10)))


if __name__ == "__main__":
    main()
