"""Exact audit for luna_one_coordinate_mutation_audit.md."""

from itertools import combinations, combinations_with_replacement
from math import gcd, isqrt, lcm
from fractions import Fraction

from check_gaussian_reflection_replacement import (
    gmul, gconj, gnorm, gcd_many, gdivexact,
)


L = 6


def pair_ok(a, b, scale=L):
    return a != b and (a * b + scale * scale) % (a - b) == 0


def divisors(n):
    n = abs(n)
    out = []
    for d in range(1, isqrt(n) + 1):
        if n % d:
            continue
        out.extend((d, -d))
        if d * d != n:
            out.extend((n // d, -n // d))
    return out


def edge_norm(q, scale=L):
    g = gcd(abs(q), scale)
    a, b = q // g, scale // g
    return (a * a + b * b) // (2 if a % 2 and b % 2 else 1)


def all_edge_norms(xs, scale=L):
    assert all(pair_ok(a, b, scale) for a, b in combinations(xs, 2))
    qs = list(xs) + [(a * b + scale * scale) // (a - b)
                     for a, b in combinations(xs, 2)]
    return [edge_norm(q, scale) for q in qs]


def conductor(xs, scale=L):
    n = 1
    for m in all_edge_norms(xs, scale):
        n = n * m // gcd(n, m)
    return n


def replacement_candidates(retained):
    y0 = retained[0]
    out = set()
    for d in divisors(y0 * y0 + L * L):
        x = y0 + d
        if x >= L and x not in retained and all(pair_ok(x, y) for y in retained):
            out.add(x)
    return sorted(out)


def check_clique(xs, scale=L):
    assert len(set(xs)) == len(xs)
    assert all(pair_ok(a, b, scale) for a, b in combinations(xs, 2))


def cmul(a, b):
    return (a[0] * b[0] - a[1] * b[1],
            a[0] * b[1] + a[1] * b[0])


def cdiv(a, b):
    den = b[0] * b[0] + b[1] * b[1]
    return ((a[0] * b[0] + a[1] * b[1]) / den,
            (a[1] * b[0] - a[0] * b[1]) / den)


def phase(x, scale):
    den = x * x + scale * scale
    return (Fraction(x * x - scale * scale, den),
            Fraction(2 * x * scale, den))


def cross(a, b):
    return a[0] * b[1] - a[1] * b[0]


def primitive_phases(qs):
    den = lcm(*(v.denominator for z in qs for v in z))
    raw = [(int(den * z[0]), int(den * z[1])) for z in qs]
    common = gcd_many(raw)
    primitive = [gdivexact(z, common) for z in raw]
    assert all(gnorm(z) == gnorm(primitive[0]) for z in primitive)
    assert gnorm(gcd_many(primitive)) == 1
    return primitive


def reflection_records(qs):
    """Exact inside test for these ordered first-quadrant source phases."""
    assert qs[0] == (1, 0)
    assert all(q[0] > 0 and q[1] >= 0 for q in qs)
    assert all(cross(a, b) > 0 for a, b in combinations(qs, 2))
    n = gnorm(primitive_phases(qs)[0])
    records = []
    for deleted in range(len(qs)):
        retained = [i for i in range(len(qs)) if i != deleted]
        for i, j in combinations_with_replacement(retained, 2):
            reflected = cdiv(cmul(qs[i], qs[j]), qs[deleted])
            new = qs[:deleted] + [reflected] + qs[deleted + 1:]
            assert len(set(new)) == len(new)
            nnew = gnorm(primitive_phases(new)[0])
            inside = cross(qs[0], reflected) >= 0 and cross(reflected, qs[-1]) >= 0
            records.append((deleted, i, j, inside, Fraction(nnew, n), reflected))
    return records


def main():
    chain = [(8, 9, 10, 12), (9, 10, 12, 18), (10, 12, 18, 27)]
    expected = [27625, 5525, 425]
    for xs, n in zip(chain, expected):
        check_clique(xs)
        assert conductor(xs) == n

    # q(27)=q(10)/q(18), checked by exact Gaussian rational arithmetic.
    def q(x):
        return phase(x, L)

    assert q(27) == cdiv(q(10), q(18))
    assert q(27) != cdiv(cmul(q(10), q(18)), q(9))

    # Include the infinity anchor q(infinity)=1. For the C1 -> C2 move, q(27)
    # lies strictly inside the old anchored span [q(infinity), q(9)]. The second
    # check uses a repeated retained index, i=j=10.
    anchor = (Fraction(1), Fraction(0))
    q27 = q(27)
    q9 = q(9)
    assert q27[1] > 0
    assert cmul((q27[0], -q27[1]), q9)[1] > 0
    q10_twice_over_q9 = cdiv(cmul(q(10), q(10)), q9)
    assert q10_twice_over_q9[1] > 0
    assert cmul((q10_twice_over_q9[0], -q10_twice_over_q9[1]), q9)[1] > 0
    assert q27 != anchor and q27 != q9 and q27 != q(10)

    # Genuine reflection: delete18, use retained infinity and10, add27.
    c1_phases = [anchor] + [q(x) for x in (18, 12, 10, 9)]
    c1_records = reflection_records(c1_phases)
    reflected_move = next(r for r in c1_records if r[:3] == (1, 0, 3))
    assert reflected_move[3] and reflected_move[4] == Fraction(1, 5)
    assert reflected_move[5] == q(27)
    assert not pair_ok(9, 27)
    assert Fraction(9 * 27 + L * L, 9 - 27) == Fraction(-31, 2)
    check_clique((18, 20, 24, 54), scale=12)
    assert conductor((18, 20, 24, 54), scale=12) == 1105
    assert [phase(x, 12) for x in (18, 20, 24, 54)] == [q(x) for x in (9, 10, 12, 27)]

    # Reconstruct the first four-point Pell witness from its Gaussian blocks.
    P, A, B, C = (-1, 2), (13, 8), (8, 5), (50, 31)
    neg = lambda z: (-z[0], -z[1])
    pell = [
        neg(gmul(gmul(gmul(P, A), B), C)),
        neg(gmul(gmul(gmul(gconj(P), A), gconj(B)), gconj(C))),
        neg(gmul(gmul(gmul(gconj(P), gconj(A)), B), gconj(C))),
        neg(gmul(gmul(gmul(gconj(P), gconj(A)), gconj(B)), C)),
    ]
    assert pell == [(16069, 10032), (16197, 9824), (16059, 10048), (16131, 9932)]
    assert gnorm(pell[0]) == 358853785
    assert gnorm(gcd_many(pell)) == 1
    ordered_pell = [pell[i] for i in (1, 3, 0, 2)]
    pell_phases = [cdiv(tuple(map(Fraction, z)), tuple(map(Fraction, ordered_pell[0])))
                   for z in ordered_pell]
    assert pell_phases == [anchor] + [phase(x, 24) for x in (7184, 3723, 3456)]
    pell_records = reflection_records(pell_phases)
    inside = [r for r in pell_records if r[3]]
    assert len(pell_records) == 24 and len(inside) == 8
    assert min(r[4] for r in inside) == 233
    assert min(r[4] for r in pell_records) == 89
    assert all(not r[3] for r in pell_records if r[4] == 89)

    final = chain[-1]
    assert replacement_candidates(list(final)) == []
    for i, x in enumerate(final):
        retained = list(final[:i] + final[i + 1:])
        for replacement in replacement_candidates(retained):
            if replacement > x and min(retained + [replacement]) > min(final):
                assert conductor(retained + [replacement]) > conductor(final)

    print("PASS: singular chain step distinguished from reflection; actual inside reflection has ratio 1/5.")
    print("PASS: Pell witness has 24 reflections, 8 inside; minima 89 overall and 233 inside.")
    print("PASS: exact clique chain, monomial identity, and terminal mutation audit")


if __name__ == "__main__":
    main()
