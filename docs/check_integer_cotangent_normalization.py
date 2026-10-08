"""Exact checks for integer cotangent gcds, radii, and reanchoring.

The proofs are in the companion notes.  These finite checks retain signs,
common integer contents, the prime two, and both Gaussian orientations.
"""

from itertools import combinations
from math import gcd, lcm

from check_least_radius_formula import (
    conj,
    exact_div,
    gcd_all,
    gcd_gaussian,
    lcm_all,
    mul,
    norm,
)


def edge_quotient(x, y, scale):
    numerator = x * y + scale * scale
    assert numerator % (x - y) == 0
    return numerator // (x - y)


def edge_norm(q, scale):
    content = gcd(q, scale)
    a, b = q // content, scale // content
    epsilon = 2 if a % 2 and b % 2 else 1
    assert (a * a + b * b) % epsilon == 0
    return (a * a + b * b) // epsilon


def primitive_tuple(xs, scale):
    fractions = []
    for x in xs:
        h = (x, scale)
        common = gcd_gaussian(h, conj(h))
        fractions.append((exact_div(h, common), exact_div(conj(h), common)))
    denominator = lcm_all([b for _, b in fractions])
    rows = [denominator]
    rows.extend(mul(exact_div(denominator, b), a) for a, b in fractions)
    assert norm(gcd_all(rows)) == 1
    assert all(norm(z) == norm(denominator) for z in rows)
    return rows


def check_clique(xs, scale):
    assert len(set(xs)) == len(xs)
    edges = list(xs)
    edges.extend(edge_quotient(x, y, scale) for x, y in combinations(xs, 2))
    rows = primitive_tuple(xs, scale)
    ordinary_radius_squared = lcm(*(edge_norm(q, scale) for q in edges))
    assert norm(rows[0]) == ordinary_radius_squared, (xs, scale)

    # All-edge data are unchanged up to sign by swapping the anchor with a.
    for a in xs:
        reflected = [a]
        reflected.extend(edge_quotient(a, x, scale) for x in xs if x != a)
        # f_a(x)=(a*x+scale^2)/(x-a), so the quotient above needs a minus.
        reflected = [a] + [-x for x in reflected[1:]]
        reflected_rows = primitive_tuple(reflected, scale)
        assert norm(reflected_rows[0]) == ordinary_radius_squared
        new_edges = list(reflected)
        new_edges.extend(
            edge_quotient(x, y, scale) for x, y in combinations(reflected, 2)
        )
        assert sorted(map(abs, edges)) == sorted(map(abs, new_edges))
    return ordinary_radius_squared


def main():
    pairs = 0
    for scale in range(1, 25):
        for x, y in combinations(range(-70, 71), 2):
            delta = x - y
            if (x * x + scale * scale) % delta:
                continue
            h, k = (x, scale), (y, scale)
            d = gcd((x * x + scale * scale) // delta, x, scale, delta)
            gaussian_gcd_norm = norm(gcd_gaussian(h, k))
            assert gaussian_gcd_norm == abs(delta) * d
            q = edge_quotient(x, y, scale)
            assert scale % d == q % d == 0
            quotient = exact_div(mul(k, conj(h)), (gaussian_gcd_norm, 0))
            sign = 1 if delta > 0 else -1
            assert quotient == (sign * q // d, sign * scale // d)
            pairs += 1

    cliques = 0
    for scale in range(1, 13):
        values = list(range(-14, 15))
        compatible = {
            (x, y): (x * y + scale * scale) % (x - y) == 0
            for x, y in combinations(values, 2)
        }
        for size in (2, 3, 4):
            for xs in combinations(values, size):
                if all(compatible[x, y] for x, y in combinations(xs, 2)):
                    check_clique(xs, scale)
                    cliques += 1

    # A sharp local example, and guaranteed larger cliques obtained by
    # clearing the finitely many pair denominators of integer parameters.
    check_clique((-18, -6, -3, 0, 2, 6, 12), 6)
    cliques += 1
    for size in range(2, 11):
        scale = lcm(*range(1, size))
        check_clique(tuple(scale * (20 + j) for j in range(size)), scale)
        cliques += 1

    u, v = 1, 0
    previous_ratio = None
    pell_checks = 0
    for exponent in range(1, 22):
        u, v = 9 * u + 20 * v, 4 * u + 9 * v
        assert u * u - 5 * v * v == 1
        if exponent % 10 != 1:
            continue
        scale = 24
        xs = (90 * u * v + 30 * v * v + 3, 96 * u * v,
              160 * v * v + 128 * u * v + 16)
        n = check_clique(xs, scale)
        aa = 1 + 10 * v * v + 2 * u * v
        bb = 1 + 10 * v * v - 2 * u * v
        cc = 13 + 130 * v * v + 38 * u * v
        assert n == 5 * aa * bb * cc
        assert min(xs) == 96 * u * v
        numerator, denominator = n * scale**4, min(xs)**4
        if previous_ratio is not None:
            old_numerator, old_denominator = previous_ratio
            assert numerator * old_denominator < old_numerator * denominator
        previous_ratio = numerator, denominator
        pell_checks += 1
        cliques += 1

    # An anchor-only ordinary lcm really can lose conjugate orientations.
    xs, scale = (4, 36), 12
    all_edges_radius_squared = check_clique(xs, scale)
    anchor_lcm = lcm(*(edge_norm(x, scale) for x in xs))
    assert anchor_lcm == 5 and all_edges_radius_squared == 25
    print(f"Exact Gaussian gcd and signed quotient checks: {pairs} pairs.")
    print(f"All-edge integer lcm and all-anchor reflection checks: {cliques} cliques.")
    print(f"Exact Pell radius comparisons: {pell_checks}.")
    print("Anchor-only counterexample: L=12, X=(4,36), lcm=5 but radius^2=25.")


if __name__ == "__main__":
    main()
