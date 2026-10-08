"""Literal full-cut fixtures for the five-row Vieta core transfer.

The invariant is constructed to vanish and has unrestricted coefficient
height. These fixtures verify identities, not small-height endpoint cases.
"""

from functools import reduce
from math import gcd

from check_five_row_gradient_cofactor_stress import (
    gaussian_gcd, mul, norm, power, split_primes,
)
from check_five_row_gradient_discriminant_core_count import (
    evaluate, independent_graphs, quadratic_coefficients,
)
from check_smaller_cut_relation_rank import add, scale


def extended_gcd(a, b):
    if b == 0:
        return abs(a), 1 if a >= 0 else -1, 0
    g, s, t = extended_gcd(b, a % b)
    return g, t, s - (a // b) * t


def bracket(p, q):
    return p[0] * q[1] - p[1] * q[0]


def check_literal_core_transfers():
    basis = independent_graphs(5, 2)
    primes = split_primes(32)
    blocks = {u: primes[u - 1] for u in range(1, 32)}
    checks = 0
    for correction_exponents in ((0, 0, 0, 0, 0), (0, 1, 2, 1, 3)):
        corrections = [power(primes[-1], e) for e in correction_exponents]
        points = []
        for i in range(5):
            p = corrections[i]
            for u in range(1, 32):
                if u & (1 << i):
                    p = mul(p, blocks[u])
            assert gcd(abs(p[0]), abs(p[1])) == 1
            assert norm(gaussian_gcd(p, (p[0], -p[1]))) == 1
            points.append(p)
        values = sum(points, ())

        # This construction enforces a numerical zero but no height bound.
        trial = {}
        for j, (_, graph) in enumerate(basis[1:]):
            trial = add(trial, scale(graph, j + 1))
        q = add(scale(trial, evaluate(basis[0][1], values)),
                scale(basis[0][1], -evaluate(trial, values)))
        content = reduce(gcd, (abs(x) for x in q.values()))
        q = {m: c // content for m, c in q.items()}
        assert q and evaluate(q, values) == 0

        for i, (x, y) in enumerate(points):
            cpoly, bpoly, apoly = quadratic_coefficients(q, i)
            other_values = values[:2 * i] + values[2 * i + 2:]
            a, b, c = [evaluate(poly, other_values)
                       for poly in (apoly, bpoly, cpoly)]
            g = reduce(gcd, (abs(a), abs(b), abs(c)))
            lam = ((2 * a * x + b * y) // (-y) if y
                   else (b * x + 2 * c * y) // x)
            assert (2 * a * x + b * y, b * x + 2 * c * y) == (-lam * y, lam * x)
            assert lam != 0

            gg, s, t = extended_gcd(x, y)
            assert gg == 1
            frame = (-t, s)
            assert bracket((x, y), frame) == 1
            frame_value = (a * frame[0] ** 2 + b * frame[0] * frame[1]
                           + c * frame[1] ** 2)
            assert gcd(abs(lam), abs(frame_value)) == g
            assert (frame_value * x - lam * frame[0]) % g == 0
            assert (frame_value * y - lam * frame[1]) % g == 0
            new = ((frame_value * x - lam * frame[0]) // g,
                   (frame_value * y - lam * frame[1]) // g)
            assert gcd(abs(new[0]), abs(new[1])) == 1
            assert a * new[0] ** 2 + b * new[0] * new[1] + c * new[1] ** 2 == 0
            assert bracket(points[i], new) != 0

            forced_content = 1
            selected = []
            for u in range(1, 32):
                inside = bool(u & (1 << i))
                other_count = bin(u & ~(1 << i)).count("1")
                exponent = max(0, 2 * other_count - (5 if inside else 4))
                forced_content *= norm(blocks[u]) ** exponent
                if ((not inside and other_count >= 2)
                        or (inside and other_count >= 3)):
                    selected.append(u)
            assert len(selected) == 16
            assert g * norm(corrections[i]) % forced_content == 0
            excess = g * norm(corrections[i]) // forced_content
            selected_core = (1, 0)
            for u in selected:
                selected_core = mul(selected_core, blocks[u])
                assert norm(gaussian_gcd(mul((excess, 0), new), blocks[u])) == norm(blocks[u])
            retained_core = gaussian_gcd(selected_core, new)
            lost_norm = norm(selected_core) // norm(retained_core)
            assert excess % lost_norm == 0

            for j in range(5):
                if j == i:
                    continue
                shared = [u for u in selected if u & (1 << j)]
                assert len(shared) == 11
                pair_core_norm = 1
                for u in shared:
                    pair_core_norm *= norm(blocks[u])
                new_gcd_norm = norm(gaussian_gcd(new, points[j]))
                assert new_gcd_norm * excess >= pair_core_norm
                q_at_j = (a * points[j][0] ** 2 + b * points[j][0] * points[j][1]
                          + c * points[j][1] ** 2)
                assert q_at_j == (g * bracket(points[i], points[j])
                                  * bracket(new, points[j]))
                assert bracket(new, points[j]) % new_gcd_norm == 0
            checks += 1
    assert checks == 10
    return checks


if __name__ == "__main__":
    checks = check_literal_core_transfers()
    print(f"PASS: {checks} literal Vieta swaps: primitive normalization and all four pair formulas.")
    print("PASS: content19 with corrections; selected16/shared11; exact lost-core norm divisibility.")
    print("The vanishing invariant has unrestricted coefficient height; no endpoint family is asserted.")
