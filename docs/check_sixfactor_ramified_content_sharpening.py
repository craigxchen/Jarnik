"""Exact bounded audit of the 2-adic content sharpening.

The proof used in the accompanying audit is the valuation case split
``delta_2 <= 1``.  This script checks that split exhaustively on a finite
residue/valuation box, and checks the resulting safe global bound
Norm(G) <= 50*P^2 on actual signed Gaussian tuples.
"""

import math
from functools import reduce

from check_joint_affine_relation_lattice import gdiv_exact, gmul, ggcd, gnorm


SUBSETS = ({4, 6}, {1, 3, 6}, {1, 4, 5}, {2, 3, 5}, {1, 2, 3, 4})


def v2(value):
    if value == 0:
        return 10**9
    value = abs(value)
    answer = 0
    while value % 2 == 0:
        value //= 2
        answer += 1
    return answer


def delta2(u, v, t):
    coefficients = (u, 2*u, v, u+v, 2*u+v, 3*u+v)
    valuations = [v2(a) for a in coefficients]
    return sum(2*min(s, t) + (s == t) - 2*s for s in valuations)


def signed_points(u, v, t):
    coefficients = (u, 2*u, v, u+v, 2*u+v, 3*u+v)
    points = []
    for subset in SUBSETS:
        z = (1, 0)
        for index, coefficient in enumerate(coefficients, 1):
            z = gmul(z, (coefficient, t if index in subset else -t))
        points.append((-z[0], -z[1]) if len(subset) % 2 else z)
    return coefficients, points


def main():
    valuation_checks = 0
    for u in range(-255, 256):
        for v in range(-255, 256):
            coefficients = (u, 2*u, v, u+v, 2*u+v, 3*u+v)
            if any(a == 0 for a in coefficients):
                continue
            for t in range(9):
                if math.gcd(math.gcd(abs(u), abs(v)), 2**t) != 1:
                    continue
                assert delta2(u, v, t) <= 1
                valuation_checks += 1

    tuple_checks = 0
    for u in range(-10, 11):
        for v in range(-20, 21):
            coefficients = (u, 2*u, v, u+v, 2*u+v, 3*u+v)
            if (any(a == 0 for a in coefficients)
                    or len(set(map(abs, coefficients))) < 6):
                continue
            P = abs(math.prod(coefficients))
            for t in range(1, 301):
                if math.gcd(math.gcd(abs(u), abs(v)), t) != 1:
                    continue
                _, points = signed_points(u, v, t)
                common = reduce(ggcd, points)
                primitive = [gdiv_exact(z, common) for z in points]
                n = gnorm(primitive[0])
                assert all(gnorm(z) == n for z in primitive)
                assert gnorm(common) <= 50*P*P

                max_chord_sq = max(
                    (primitive[i][0] - primitive[j][0])**2
                    + (primitive[i][1] - primitive[j][1])**2
                    for i in range(5) for j in range(i)
                )
                # Exact square of the normalized-chord inequality
                # C_chord >= 2*sqrt(P/|G|).
                assert max_chord_sq**2 * gnorm(common) >= 16*P*P*n
                tuple_checks += 1

    print(f"Checked {valuation_checks} exact 2-adic cases and "
          f"{tuple_checks} signed tuples; delta_2 <= 1 and "
          "Norm(G) <= 50 P^2 verified.")


if __name__ == "__main__":
    main()
