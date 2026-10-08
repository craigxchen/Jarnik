"""Exact polynomial and arithmetic checks of the five-point cubic family."""

from functools import reduce
from itertools import combinations
from math import gcd

from check_complement_reflection_invariant_relations import (
    ONE, conjugate, exact_div, gaussian_gcd, mul as gmul,
    norm, product,
)
from check_scaled_split_contact_profile_q2 import add, mul, scale, trim


SUBSETS = [{4, 6}, {1, 3, 6}, {1, 4, 5}, {2, 3, 5}, {1, 2, 3, 4}]
PQ = [(27, 108), (56, 1008), (144, 3600), (50, 900), (225, 3600),
      (200, 3600), (21, 0), (108, 108), (112, 1008), (25, 900)]


def evaluate(p, t):
    result = 0
    for c in reversed(p):
        result = result * t + c
    return result


def forward_point(subset):
    x, y = [(-1)**len(subset)], [0]
    for a in range(1, 7):
        sign = 1 if a in subset else -1
        x, y = (add(mul(x, [0, a]), scale(y, -sign)),
                add(mul(y, [0, a]), scale(x, sign)))
    return x, y


def reverse_degree_six(p):
    return trim(list(reversed(p + [0] * (7-len(p)))))


def triangle_polynomial(points, triple):
    i, j, k = triple
    x, y = points[i]
    xx, yy = points[j]
    xxx, yyy = points[k]
    return add(mul(add(xx, scale(x, -1)), add(yyy, scale(y, -1))),
               scale(mul(add(yy, scale(y, -1)), add(xxx, scale(x, -1))), -1))


def triangle(points, triple):
    i, j, k = triple
    a = (points[j][0]-points[i][0], points[j][1]-points[i][1])
    b = (points[k][0]-points[i][0], points[k][1]-points[i][1])
    return a[0]*b[1]-a[1]*b[0]


def prime_factors(n):
    out = set()
    d = 2
    while d*d <= n:
        while n % d == 0:
            out.add(d)
            n //= d
        d += 1
    if n > 1:
        out.add(n)
    return out


def main():
    assert all(sum(s) == 10 for s in SUBSETS)
    assert all(any(a in s for s in SUBSETS) and any(a not in s for s in SUBSETS)
               for a in range(1, 7))
    coeffs = [60//a for a in range(1, 7)]
    possible_primes = {2}
    for a, b in combinations(coeffs, 2):
        possible_primes |= prime_factors(a+b) | prime_factors(abs(a-b))
    assert possible_primes == {2, 3, 5, 7, 11}
    assert {p for p in possible_primes if p == 2 or p % 4 == 1} == {2, 5}

    forward = [forward_point(s) for s in SUBSETS]
    reverse = [(reverse_degree_six(x), reverse_degree_six(y)) for x, y in forward]
    triples = list(combinations(range(5), 3))
    target_norm = [1]
    for a in range(1, 7):
        target_norm = mul(target_norm, [a*a, 0, 1])
    for x, y in reverse:
        assert add(mul(x, x), mul(y, y)) == target_norm
    for triple, (p, q) in zip(triples, PQ):
        assert triangle_polynomial(forward, triple) == trim([0]*9+[-960*p, 0, -960*q])
        assert triangle_polynomial(reverse, triple) == trim([0, -960*q, 0, -960*p])
    assert triangle_polynomial(forward, (1, 2, 3)) == [0]*9+[-20160]

    # For T=6s the q-polynomial gcd is exactly 36 for every integer s:
    # the divided entries satisfy 4*q012-q124=9 and q024=1 mod3.
    for s in range(1, 121):
        tt = 6*s
        assert reduce(gcd, (p*tt*tt+q for p, q in PQ)) == 36

    for k in list(range(1, 41)) + [100, 1000]:
        tt = 600*k
        raw = [(evaluate(x, tt), evaluate(y, tt)) for x, y in reverse]
        by_product = [gmul(((-1)**len(s), 0), product((a, tt if a in s else -tt)
                       for a in range(1, 7))) for s in SUBSETS]
        assert raw == by_product
        points = [exact_div(z, (720, 0)) for z in raw]
        assert len(set(points)) == 5
        assert norm(reduce(gaussian_gcd, points)) == 1
        nn = norm(points[0])
        assert all(norm(z) == nn for z in points)
        determinants = [triangle(points, triple) for triple in triples]
        assert all(d < 0 for d in determinants)
        gg = abs(reduce(gcd, determinants))
        assert gg == tt//15
        hh = max(map(abs, determinants))//gg
        assert hh == 25*tt*tt//4+100
        # The affine-normalized conic and actual Gram scale.
        x = tt*tt
        aa = 4*(x+25)*(x+36)
        bb = 10*(x*x+43*x+360)
        cc = 25*(x+9)*(x+16)
        assert aa*cc-bb*bb == (180*tt)**2
        assert aa+cc-2*bb == 9*x*(x+1)
        basis = [(points[j][0]-points[0][0], points[j][1]-points[0][1])
                 for j in (1, 2)]
        assert norm(basis[0])*3600 == (x+4)*aa
        assert norm(basis[1])*3600 == (x+4)*cc
        assert sum(a*b for a, b in zip(*basis))*3600 == (x+4)*bb
        for un, vn, den in [(-50*(x+18), 56*(x+18), 27*(x+4)),
                            (-25*(x+16), 16*(x+25), 3*(x+4))]:
            assert aa*un*un+2*bb*un*vn+cc*vn*vn-aa*un*den-cc*vn*den == 0
        assert all(z[0] < 0 and z[1] < 0 for z in points)
        assert all(points[j][0] > points[j+1][0] and points[j][1] < points[j+1][1]
                   for j in range(4))
        chord = (points[4][0]-points[0][0], points[4][1]-points[0][1])
        dd = norm(chord)
        assert dd*36 == tt*tt*(tt**4+41*tt*tt+400)
        # d/(2R)<=1/1000 and d/sqrt(R)<=2sqrt(5)*sqrt(1001/1000).
        assert dd*1000**2 <= 4*nn
        assert dd**2 * 1000**2 <= 20**2 * 1001**2 * nn
    # The two elementary bounds imply C_arc<5 by arcsin(u)<=u/sqrt(1-u²).
    assert 20*1001*1000 < 25*(1000**2-1)
    print('PASS: all forward and reversed polynomial norm/triangle identities.')
    print('PASS: exact t^9 divisor and exhaustive possible-prime content certificate.')
    print('PASS: 42 primitive endpoint specializations, exact joint height, and C_arc<5.')
    print('Formulas give R/H_shape² -> infinity and C_arc -> 2sqrt(5).')


if __name__ == '__main__':
    main()
