"""Exact Gaussian cost of adjoining a reflected circle configuration."""

from itertools import combinations
from math import gcd, lcm

from check_integer_cotangent_normalization import primitive_tuple
from check_least_radius_formula import (
    conj, exact_div, gcd_all, gcd_gaussian, mul, norm,
)


def audit(rows, i, j):
    n = norm(rows[0])
    assert norm(gcd_all(rows)) == 1
    assert all(norm(z) == n for z in rows)
    z0, za = rows[i], rows[j]
    common = gcd_gaussian(za, conj(z0))
    alpha = exact_div(za, common)
    beta = exact_div(conj(z0), common)
    assert norm(gcd_gaussian(alpha, beta)) == 1
    d = norm(beta)
    assert norm(alpha) == d
    union = [mul(beta, z) for z in rows]+[mul(alpha, conj(z)) for z in rows]
    assert norm(gcd_all(union)) == 1
    assert all(norm(z) == n*d for z in union)
    edge_common = gcd_gaussian(za, z0)
    edge_norm = norm(exact_div(z0, edge_common))
    contents = gcd(gcd(*z0), gcd(*za))
    assert d*edge_norm*contents**2 == n
    chord = (za[0]-z0[0], za[1]-z0[1])
    g = gcd(*chord)
    u, v = chord[0]//g, chord[1]//g
    eps = 2 if u % 2 and v % 2 else 1
    assert d*eps == u*u+v*v
    assert norm(gcd_gaussian(alpha, conj(beta))) == d
    return d, len(set(union))


def main():
    cases = 0
    for size in range(2, 9):
        scale = lcm(*range(1, size))
        for offset in (0, 4, 20):
            rows = primitive_tuple([scale*(offset+j) for j in range(size)], scale)
            for i, j in combinations(range(len(rows)), 2):
                audit(rows, i, j)
                cases += 1
    u, v = 1, 0
    pell = 0
    for exponent in range(1, 22):
        u, v = 9*u+20*v, 4*u+9*v
        if exponent % 10 != 1:
            continue
        xs = [90*u*v+30*v*v+3, 96*u*v, 160*v*v+128*u*v+16]
        rows = primitive_tuple(xs, 24)
        d, count = audit(rows, 0, 2)
        a = 1+10*v*v+2*u*v
        b = 1+10*v*v-2*u*v
        cp = 13+130*v*v+38*u*v
        assert min(xs) == xs[1]
        assert norm(rows[0]) == 5*a*b*cp
        assert a*b == 16*u*u*v*v+1
        assert gcd(*rows[0]) == gcd(*rows[2]) == 1
        assert d == 5*cp and count == 6
        pell += 1
    print(f'PASS: {cases} arbitrary tuple/pair mirror costs; {pell} exact Pell unions.')


if __name__ == '__main__':
    main()
