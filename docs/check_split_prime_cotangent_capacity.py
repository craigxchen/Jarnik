"""Exact positive sharpness examples for the split-prime capacity bound."""

from itertools import combinations
from math import gcd, lcm


def valuation(value, prime):
    assert value
    power = 0
    while value % prime == 0:
        value //= prime
        power += 1
    return power


def norm_lcm(coordinates, denominator):
    edges = list(coordinates)
    for x, y in combinations(coordinates, 2):
        quotient, remainder = divmod(x * y + denominator ** 2, y - x)
        assert remainder == 0
        edges.append(quotient)
    norms = []
    for edge in edges:
        common = gcd(edge, denominator)
        a, b = edge // common, denominator // common
        epsilon = 2 if a % 2 and b % 2 else 1
        norms.append((a * a + b * b) // epsilon)
    return lcm(*norms), norms


def main():
    signed = [-18, -6, -3, 0, 2, 6, 12]
    modulus = lcm(*(y - x for x, y in combinations(signed, 2)))
    assert modulus == 360 and lcm(modulus, 25) == 1800
    radius, norms = norm_lcm(signed, 6)
    assert radius == 5 and set(norms) == {1, 5}
    cases = [(6, signed, 5),
             (12, [264, 288, 294, 300, 304, 312, 324], 2409083501772215645),
             (6, [x + 1800 for x in signed],
              854184548609224587096599042100005)]
    for denominator, coordinates, expected in cases:
        radius, _ = norm_lcm(coordinates, denominator)
        assert radius == expected
        assert valuation(radius, 5) == 1 and valuation(denominator, 5) == 0
        assert len(coordinates) + 1 == (valuation(radius, 5) + 1) * 4
        for prime in (5, 13, 17, 29, 37):
            assert len(coordinates) + 1 <= ((valuation(radius, prime) + 1) *
                    (prime - 1) * prime ** valuation(denominator, prime))
    for multiplier in (1, 2, 3, 5, 13, 25, 101):
        coordinates = [x + 1800 * multiplier for x in signed]
        assert min(coordinates) >= 6
        radius, _ = norm_lcm(coordinates, 6)
        assert valuation(radius, 5) == 1
    print("PASS: all three exact norm lcms and sharp split-prime capacities")
    print("PASS: seven positive fixed-denominator translations, retaining v_5(N)=1")


if __name__ == "__main__":
    main()
