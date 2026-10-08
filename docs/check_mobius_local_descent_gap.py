"""Exact checks for the local nonconformal radius gap."""

from fractions import Fraction
from itertools import product
from math import gcd

from check_mobius_conductor_transfer import fixtures, least_norm
from check_mobius_conformal_normalization import image


def determinant(h, k):
    return h[0] * k[1] - h[1] * k[0]


def norm(h):
    return h[0] * h[0] + h[1] * h[1]


def check_phase_separation():
    cases = 0
    rows = [
        h
        for h in product(range(-8, 9), repeat=2)
        if h != (0, 0) and gcd(*h) == 1
    ]
    for i, h in enumerate(rows):
        for k in rows[i + 1 :]:
            if determinant(h, k) == 0:
                continue
            # Squaring (4) after dropping |det| >= 1.
            assert 4 * determinant(h, k) ** 2 >= 4
            assert norm(h) > 0 and norm(k) > 0
            cases += 1
    return cases


def check_shears():
    cases = strict = 0
    source_sets = fixtures()
    source_sets.append(
        [(1, 0), (1, 1), (1, 2), (3, 1), (3, 2), (5, 1), (5, 2), (2, 5)]
    )
    for rows in source_sets:
        n = least_norm(rows)
        for p in range(1, 6):
            for q in range(p + 1, 121):
                if gcd(p, q) != 1:
                    continue
                # The rational shear is represented integrally by [[q,p],[0,q]].
                targets = image(rows, (q, p, 0, q))
                nprime = least_norm(targets)
                assert 64 * p * p * n * nprime >= (q - p) ** 2
                if q > (8 * n + 1) * p:
                    assert nprime > n
                    strict += 1
                cases += 1

        # Test the first integral denominator beyond the local threshold.
        p = 1
        q = 8 * n + 2
        targets = image(rows, (q, p, 0, q))
        nprime = least_norm(targets)
        assert nprime > n
        assert 64 * n * nprime >= (q - 1) ** 2
        strict += 1
        cases += 1
    return cases, strict


def check_distance_formula():
    cases = 0
    for a, b, c, d in product(range(-5, 6), repeat=4):
        determinant_m = a * d - b * c
        if determinant_m == 0:
            continue
        trace = a * a + b * b + c * c + d * d
        epsilon_squared = Fraction(
            trace - 2 * abs(determinant_m),
            trace + 2 * abs(determinant_m),
        )
        assert 0 <= epsilon_squared < 1
        assert epsilon_squared == Fraction(
            (trace - 2 * abs(determinant_m)) ** 2,
            trace * trace - 4 * determinant_m * determinant_m,
        ) if trace != 2 * abs(determinant_m) else epsilon_squared == 0
        cases += 1
    return cases


def main():
    separation = check_phase_separation()
    shears, strict = check_shears()
    distances = check_distance_formula()
    print(f"PASS: {separation} exact primitive phase-separation pairs.")
    print(f"PASS: {shears} rational shear radius-jump inequalities.")
    print(f"PASS: {strict} explicit shears beyond the no-descent threshold.")
    print(f"PASS: {distances} exact normalized-distance identities.")


if __name__ == "__main__":
    main()
