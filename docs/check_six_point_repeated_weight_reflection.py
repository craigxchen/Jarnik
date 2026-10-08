"""Exact repeated-weight reflection and single-opposite-pair checks."""

from fractions import Fraction as F
from functools import reduce
from itertools import combinations_with_replacement

from check_integer_cotangent_spectral import (
    ONE, add, conj, div, mul, neg, power, primitive_vector,
)
from check_six_point_fixed_moment_elliptic import P, ec_add, rational_sqrt


def moments(rows, weights):
    for exponent in range(3):
        assert reduce(add, (mul((u, 0), power(z, exponent))
                            for u, z in zip(weights, rows)), (0, 0)) == (0, 0)


def main():
    point = None
    checks = 0
    for _ in range(15):
        point = ec_add(point, P)
        x, y = point
        q, p = rational_sqrt(x), y/(99-x)
        first = (1-p*q, p+q)
        second = (1+p*q, p-q)
        third = (1-q*q, 10*p)
        raw = [first, conj(first), second, conj(second), third, conj(third)]
        rows = [div(z, raw[0]) for z in raw]
        for permutation in ((0,1,2,3,4,5), (0,3,2,1,4,5)):
            zs = [rows[i] for i in permutation]
            assert len(set(zs)) == 6
            assert all(mul(z, conj(z)) == ONE for z in zs)
            moments(zs, (5,-5,5,-5,-1,1))
            t = mul(zs[4], zs[5])
            matchings = [mul(zs[0], zs[j]) == t and mul(zs[2], zs[k]) == t
                         for j, k in ((1,3), (3,1))]
            assert sum(matchings) == 1
            phi = [add(z, div(t, z)) for z in zs]
            assert add(phi[0],phi[2]) == add(phi[1],phi[3])
            assert add(power(phi[0],2),power(phi[2],2)) == add(power(phi[1],2),power(phi[3],2))
            checks += 1

    nodes = [ONE]+[div((x,1),(x,-1)) for x in (2,3,5,12,17)]
    raw = []
    for i, z in enumerate(nodes):
        denominator = reduce(mul, (add(z,neg(w)) for j,w in enumerate(nodes) if i!=j), ONE)
        raw.append(div(power(z,2),denominator))
    central, _ = primitive_vector(raw)
    expected = (126,-7,50,-169,841,-841)
    unit = div(central[0],(expected[0],0))
    assert unit in (ONE,(-1,0),(0,1),(0,-1))
    assert central == [mul(unit,(u,0)) for u in expected]
    moments(nodes, expected)
    products = [mul(nodes[i],nodes[j]) for i,j in combinations_with_replacement(range(6),2)]
    assert len(set(products)) == 21
    print(f'PASS: {checks} rational reflection configurations, both matchings, '
          'and a collision-free single-opposite-pair example.')


if __name__ == '__main__':
    main()
