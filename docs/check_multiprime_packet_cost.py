#!/usr/bin/env python3
"""Finite exact checks for multiprime packet splitting and its cost limits.

These checks supplement the proof in multiprime_packet_cost_obstruction.md;
they do not prove a uniform circle bound or an unbounded cost asymptotic.
"""

from fractions import Fraction
from itertools import combinations
from math import gcd, prod
from random import Random

from check_mixed_affine_prime_amplification import (
    divisors, euler_phi, mobius, ramanujan, row_widths,
)


def moment(row, k):
    return sum(b * ramanujan(e, k) for e, b in row.items())


def split_row(row, primes):
    """The fixed triangular linear map, plus unchanged deep packets."""
    packets = {p: {} for p in primes}
    cubes = {}
    for e, b in row.items():
        deep = next((p for p in primes if e % (p * p) == 0), None)
        if deep is not None:
            target = packets[deep]
            target[e // deep] = target.get(e // deep, 0) + b
            continue
        mask = 0
        s = e
        for j, p in enumerate(primes):
            if s % p == 0:
                mask |= 1 << j
                s //= p
        cubes.setdefault(s, {})[mask] = b
    alternating = {}
    for s, cube in cubes.items():
        alternating[s] = sum((-1) ** bin(m).count("1") * b for m, b in cube.items())
        for j, p in enumerate(primes):
            for tail in range(1 << (len(primes) - j - 1)):
                tail_mask = tail << (j + 1)
                coefficient = sum(
                    (-1) ** bin(prefix).count("1")
                    * cube.get(prefix | (1 << j) | tail_mask, 0)
                    for prefix in range(1 << j)
                )
                index = s * prod(
                    primes[l] for l in range(j + 1, len(primes))
                    if tail_mask & (1 << l)
                )
                target = packets[p]
                target[index] = target.get(index, 0) + coefficient
    return packets, alternating


def normalized(rows):
    indices = set().union(*(row.keys() for row in rows))
    minima = {e: min(row.get(e, 0) for row in rows) for e in indices}
    return [{e: row.get(e, 0) - minima[e] for e in indices} for row in rows]


def family_cost(rows):
    return sum(euler_phi(e) * w for e, w in row_widths(rows).items())


def check_split_families():
    rng = Random(527)
    families = identities = low_jets = 0
    for count in range(1, 5):
        primes = (3, 5, 7, 11)[:count]
        for _ in range(40):
            rows = []
            for _ in range(5):
                row = {}
                for s in (1, 13):
                    cube = {mask: rng.randrange(-3, 4) for mask in range(1, 1 << count)}
                    cube[0] = -sum((-1) ** bin(mask).count("1") * b for mask, b in cube.items())
                    for mask, b in cube.items():
                        e = s * prod(p for j, p in enumerate(primes) if mask & (1 << j))
                        row[e] = b
                    for p in primes:
                        row[s * p * p] = rng.randrange(-3, 4)
                if count > 1:
                    row[primes[0] ** 2 * primes[1] ** 2] = rng.randrange(-3, 4)
                # A nonzero discarded part with all odd moments below 13 zero.
                u = rng.randrange(-3, 4)
                row[1] += u
                row[13] += u
                rows.append(row)
            rows = normalized(rows)
            results = [split_row(row, primes) for row in rows]
            new_cost = sum(p * family_cost([r[0][p] for r in results]) for p in primes)
            assert new_cost <= 3 * family_cost(rows)
            families += 1
            for i, j in combinations(range(len(rows)), 2):
                delta = {e: rows[i].get(e, 0) - rows[j].get(e, 0) for e in rows[i]}
                for k in range(1, 36, 2):
                    output = sum(
                        p * (moment(results[i][0][p], k // p) - moment(results[j][0][p], k // p))
                        for p in primes if k % p == 0
                    )
                    discarded = moment(results[i][1], k) - moment(results[j][1], k)
                    assert moment(delta, k) - output == discarded
                    identities += 1
                    if k < 13:
                        assert discarded == 0
                        assert moment(delta, k) == output
                        if gcd(k, prod(primes)) == 1:
                            assert moment(delta, k) == 0
                        low_jets += 1
    return families, identities, low_jets


def check_gauge_obstruction():
    checks = 0
    for p, q in combinations((3, 5, 7, 11, 13), 2):
        claimed = p * q - abs(p - q)
        old = 1 + (p - 1) * (q - 1)
        assert claimed == old + 2 * (min(p, q) - 1)
        for scale in (1, 2, 3):
            costs = []
            for numerator in range(-16, 33):
                g = Fraction(numerator, 8) * scale
                cost = (p * q + p - q) * abs(scale - g) + (p * q - p + q) * abs(g)
                assert cost >= scale * claimed
                costs.append(cost)
                checks += 1
            assert min(costs) == scale * claimed
    assert 2 * (1 + euler_phi(35)) == 50
    assert 2 * (35 - abs(5 - 7)) == 66
    return checks


def check_base_expansion():
    checks = 0
    ratios = []
    for count in range(1, 9):
        primes = (3, 5, 7, 11, 13, 17, 19, 23)[:count]
        P = prod(primes)
        ds = divisors(P)
        base_coefficients = {d: mobius(P // d) for d in ds if d > 1}
        for e in ds:
            exponent = sum(v for d, v in base_coefficients.items() if d % e == 0)
            expected = (1 if e == P else 0) - (mobius(P) if e == 1 else 0)
            assert exponent == expected
            checks += 1
        cost = sum(d * abs(v) for d, v in base_coefficients.items())
        assert cost == prod(p + 1 for p in primes) - 1
        ratios.append(Fraction(cost, euler_phi(P) + 1))
    # This finite growth check does not replace the Euler-product proof.
    assert all(a < b for a, b in zip(ratios, ratios[1:]))
    return checks, ratios[-1]


def main():
    families, identities, low_jets = check_split_families()
    gauge_checks = check_gauge_obstruction()
    expansion_checks, last_ratio = check_base_expansion()
    print(f"PASS: {families} five-row families; {identities} exact discarded-moment identities; {low_jets} preserved low jets")
    print(f"PASS: {gauge_checks} exact gauge checks, including the parity-valid cost increase 50 to 66")
    print(f"PASS: {expansion_checks} divisor-exponent identities; eight-prime base-refinement ratio {last_ratio}")


if __name__ == '__main__':
    main()
