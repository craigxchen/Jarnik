"""Symbolic linear-leading audit and exact primitive Ptolemy comparisons."""

from fractions import Fraction
from itertools import combinations
from math import gcd

from check_intrinsic_cotangent_scale_triangle_divisibility import primitive_phases
from check_least_radius_formula import exact_div, gcd_all, mul, norm


def sub(a, b):
    return a[0]-b[0], a[1]-b[1]


def determinant(points, i, j, k):
    a, b = sub(points[j], points[i]), sub(points[k], points[i])
    return a[0]*b[1]-a[1]*b[0]


def add_poly(a, b, scalar=1):
    out = dict(a)
    for monomial, value in b.items():
        out[monomial] = out.get(monomial, 0)+scalar*value
        if not out[monomial]:
            del out[monomial]
    return out


def mul_poly(a, b):
    out = {}
    for ma, va in a.items():
        for mb, vb in b.items():
            key = tuple(x+y for x, y in zip(ma, mb))
            out[key] = out.get(key, 0)+va*vb
    return {m: v for m, v in out.items() if v}


def symbolic():
    count = 0
    for m in range(3, 10):
        variables = [{}]
        for i in range(m-1):
            exponent = tuple(int(j == i) for j in range(m-1))
            variables.append({exponent: 1})

        def cubic(i, j, k):
            return mul_poly(mul_poly(add_poly(variables[j], variables[i], -1),
                                     add_poly(variables[k], variables[i], -1)),
                            add_poly(variables[k], variables[j], -1))

        anchored = {(i, j): cubic(0, i, j)
                    for i, j in combinations(range(1, m), 2)}
        for i, j, k in combinations(range(1, m), 3):
            reduced = add_poly(add_poly(anchored[i, j], anchored[i, k], -1),
                               anchored[j, k])
            assert cubic(i, j, k) == reduced
            count += 1
        for i, j in anchored:
            monomial = tuple(1 if k == i-1 else 2 if k == j-1 else 0
                             for k in range(m-1))
            for pair, polynomial in anchored.items():
                assert polynomial.get(monomial, 0) == int(pair == (i, j))
                count += 1
    return count


def primitive_ptolemy(points):
    E = {(i, j): determinant(points, 0, i, j)
         for i, j in combinations(range(1, 5), 2)}
    assert all(value > 0 for value in E.values())
    products = (E[1, 2]*E[3, 4], E[1, 4]*E[2, 3], E[1, 3]*E[2, 4])
    assert products[0]+products[1] == products[2]
    content = gcd(gcd(products[0], products[1]), products[2])
    primitive = tuple(value//content for value in products)
    chord_products = [mul(sub(points[j], points[i]), sub(points[l], points[k]))
                      for i, j, k, l in ((1, 2, 3, 4), (1, 4, 2, 3), (1, 3, 2, 4))]
    assert tuple(a+b for a, b in zip(chord_products[0], chord_products[1])) == chord_products[2]
    common = gcd_all(chord_products)
    quotients = [exact_div(z, common) for z in chord_products]
    matches = []
    for unit in ((1, 0), (-1, 0), (0, 1), (0, -1)):
        adjusted = [mul(unit, z) for z in quotients]
        if all(z[1] == 0 and z[0] > 0 for z in adjusted):
            matches.append(tuple(z[0] for z in adjusted))
    assert matches == [primitive]


def actual_checks():
    count = 0
    for n in range(2, 15):
        points = primitive_phases([(n, j) for j in range(5)])
        N = norm(points[0])
        e = Fraction(1, n)
        actual = Fraction(determinant(points, 0, 1, 2)-determinant(points, 1, 2, 3), N)
        expected = 72*e**5/((1+e*e)*(1+4*e*e)*(1+9*e*e))
        assert actual == expected > 0
        for i, j in combinations(range(1, 5), 2):
            expected_E = 4*e**3*i*j*(j-i)/((1+e*e*i*i)*(1+e*e*j*j))
            assert Fraction(determinant(points, 0, i, j), N) == expected_E
        primitive_ptolemy(points)
        count += 1
    fixtures = [
        [(-110, -35), (-109, -38), (-98, -61), (-94, -67), (-86, -77)],
        [(4, -33), (9, -32), (12, -31), (23, -24), (24, -23)],
    ]
    for points in fixtures:
        primitive_ptolemy(points)
        count += 1
    return count


if __name__ == '__main__':
    print('symbolic reduction and unique-monomial checks:', symbolic())
    print('fully primitive circle/Ptolemy fixtures:', actual_checks())
