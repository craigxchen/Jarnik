"""Exact all-row pair-selected obstruction on a bounded-endpoint family."""

from fractions import Fraction as Q
from itertools import combinations
from math import gcd, lcm, prod

from check_four_row_rational_quartic import G, trim, multiply
from check_endpoint_swap_content_and_overlap import (
    primitive_radius, all_edge_radius,
)


def divide_scalar(a, b):
    b = G(b)
    return a*b.conjugate()*G(1/(b.a*b.a+b.b*b.b))


def add(a, b):
    return trim([(a[k] if k < len(a) else G())
                 +(b[k] if k < len(b) else G())
                 for k in range(max(len(a), len(b)))])


def scale(a, c):
    return trim([v*c for v in a])


def conjugate(a):
    return [v.conjugate() for v in a]


def divrem(a, b):
    a = a[:]
    out = [G()]*max(1, len(a)-len(b)+1)
    while len(a) >= len(b) and a != [G()]:
        shift = len(a)-len(b)
        coefficient = divide_scalar(a[-1], b[-1])
        out[shift] = coefficient
        a = add(a, [G()]*shift+scale(b, -coefficient))
    return trim(out), a


def exact_divide(a, b):
    out, remainder = divrem(a, b)
    assert remainder == [G()]
    return out


def monic(a):
    return scale(a, divide_scalar(G(1), a[-1]))


def polynomial_gcd(a, b):
    while b != [G()]:
        a, b = b, divrem(a, b)[1]
    return monic(a)


def primitive(a):
    return monic(exact_divide(a, polynomial_gcd(a, conjugate(a))))


def polynomial_lcm(a, b):
    return monic(multiply(exact_divide(a, polynomial_gcd(a, b)), b))


def denominator_lcm(rows):
    answer = [G(1)]
    for row in rows:
        answer = polynomial_lcm(answer, primitive(row))
    return answer


def evaluate(a, t):
    answer = G()
    for coefficient in reversed(a):
        answer = answer*t+coefficient
    return answer


def integer_row(a, t):
    z = evaluate(a, t)
    denominator = lcm(z.a.denominator, z.b.denominator)
    u, v = int(denominator*z.a), int(denominator*z.b)
    content = gcd(u, v)
    return u//content, v//content


def check_circle_tuple(rows, denominator):
    points = [multiply(exact_divide(conjugate(denominator), conjugate(row)),
                       row) for row in rows]
    common = points[0]
    for point in points[1:]:
        common = polynomial_gcd(common, point)
    assert common == [G(1)]
    common_norm = multiply(denominator, conjugate(denominator))
    assert all(multiply(z, conjugate(z)) == common_norm for z in points)


def main():
    cs = [404, 428, 464, 484, 524]
    es = [720, -720, -720, -720, 720]
    ds = [85, 61, 25, 5, -35]
    fs = [1164, 396, 324, -276, -1236]
    source_points = [[G(e), G(0, -f), G(-c), G(0, -d),
                      G(-45), G(0, -1), G(-1)]
                     for c, e, d, f in zip(cs, es, ds, fs)]
    radius_norm = multiply(source_points[0], conjugate(source_points[0]))
    assert all(multiply(z, conjugate(z)) == radius_norm for z in source_points)
    rows = [[G(1)]]
    for z in source_points[1:]:
        relative = multiply(z, conjugate(source_points[0]))
        rows.append(primitive(add(radius_norm, conjugate(relative))))
    expected = [
        [G(1)],
        [G(0, 12), G(13), G(), G(1)],
        [G(0, 30), G(31), G(), G(1)],
        [G(0, 720), G(444), G(0, 40), G(45), G(), G(1)],
        [G(-36), G(0, 60), G(25), G(), G(1)],
    ]
    assert rows == expected
    B = rows[4]
    real = [G(c.a) for c in B]
    imag = [G(c.b) for c in B]
    us = [[G(c.a) for c in row] for row in rows]
    vs = [[G(c.b) for c in row] for row in rows]
    transverse = [add(multiply(imag, u), scale(multiply(real, v), -1))
                  for u, v in zip(us, vs)]
    expected_d = [
        [432, 0, 480, 0, 48],
        [1080, 0, 1110, 0, 30],
        [25920, 0, 10080, 0, 980, 0, 20],
    ]
    assert transverse[1:4] == [[G(c) for c in d] for d in expected_d]
    source_lcm = denominator_lcm(rows)
    assert source_lcm == [G(-720), G(0, 1164), G(404), G(0, 85),
                          G(45), G(0, 1), G(1)]
    check_circle_tuple(rows, source_lcm)
    # The fixed-parameter classification uses only 1, B, C_3.
    expected_B = [G(1)]
    for a in (1, 2, 3, -6):
        expected_B = multiply(expected_B, [G(0, a), G(1)])
    assert B == expected_B
    d3_factor = [G(20)]
    for a in (4, 9, 36):
        d3_factor = multiply(d3_factor, [G(a), G(), G(1)])
    assert transverse[3] == d3_factor
    assert trim(vs[3]) == [G(720), G(), G(40)]
    assert divrem(transverse[3], trim(vs[3]))[1] == [G(45360)]
    fixed_common = [G(1)]
    for a in (2, 3, -6):
        fixed_common = multiply(fixed_common, [G(0, a), G(1)])
    assert polynomial_gcd(B, transverse[3]) == fixed_common
    coefficient = multiply(vs[3], B)
    assert coefficient[0] == G(-25920)
    assert transverse[3][0] == G(25920)
    assert add(coefficient, scale(conjugate(coefficient), -1)) == \
        multiply([G(), G(0, 120)], vs[3])
    assert transverse[3][-1] == G(20) and coefficient[-1] == G(40)
    for alignment in (Q(1, 2), Q(2), Q(4), Q(3, 7)):
        C3 = add(transverse[3], scale(coefficient, alignment))
        assert polynomial_gcd(C3, conjugate(C3)) == [G(1)]
        assert polynomial_gcd(B, C3) == fixed_common
        fixed_subset = [rows[0], B, primitive(C3)]
        fixed_lcm = denominator_lcm(fixed_subset)
        assert len(fixed_lcm)-1 == 7
        check_circle_tuple(fixed_subset, fixed_lcm)
    data = [
        ((1, 2),
         [G(0, Q(80, 3)), G(40), G(0, Q(20, 3)), G(21), G(), G(1)],
         [G(0, 1), G(1)],
         [G(Q(80, 3)), G(0, -Q(40, 3)), G(20), G(0, -1), G(1)], 10),
        ((1, 3),
         [G(0, 360), G(240), G(0, 20), G(29), G(), G(1)],
         [G(-6), G(0, 5), G(1)],
         [G(0, -60), G(10), G(0, -5), G(1)], 9),
        ((2, 3),
         [G(0, 960), G(544), G(0, Q(160, 3)), G(Q(140, 3)), G(), G(1)],
         [G(12), G(0, -4), G(1)],
         [G(0, 80), G(Q(56, 3)), G(0, 4), G(1)], 9),
    ]
    targets = []
    for (i, j), K, common, residual, degree in data:
        transformed = []
        for h in (1, 2, 3):
            F = add(multiply(multiply(transverse[i], transverse[j]), vs[h]),
                    multiply(multiply(multiply(vs[i], vs[j]), transverse[h]), B))
            transformed.append(primitive(F))
        assert transformed[i-1] == rows[j]
        assert transformed[j-1] == rows[i]
        h = next(h for h in (1, 2, 3) if h not in (i, j))
        assert transformed[h-1] == K
        assert polynomial_gcd(K, conjugate(K)) == [G(1)]
        subset = denominator_lcm([rows[0], B, rows[i], rows[j]])
        assert subset == source_lcm
        assert polynomial_gcd(subset, K) == common
        assert exact_divide(K, common) == residual
        target_rows = [rows[0], B]+transformed
        target_lcm = denominator_lcm(target_rows)
        assert target_lcm == multiply(source_lcm, residual)
        assert len(target_lcm)-1 == degree
        check_circle_tuple(target_rows, target_lcm)
        targets.append(target_rows)
    cases = 0
    for t in (600, 1200, 2400, 4800):
        source = [integer_row(row, t) for row in rows]
        source_n = primitive_radius(source)
        assert source_n == prod(t*t+a*a for a in range(1, 7))//720**2
        assert source_n == all_edge_radius(source)
        P, Qvalue = source[4]
        assert all(Qvalue*u-P*v > 0 and v > 0 for u, v in source[1:4])
        for target_rows in targets:
            target = [integer_row(row, t) for row in target_rows]
            target_n = primitive_radius(target)
            assert target_n == all_edge_radius(target)
            assert target_n > source_n*t**4
            cases += 1
    print('Passed: original five-point phases; three exact pair-selected maps; '
          'full Gaussian polynomial lcms of degrees 10, 9, 9; primitive '
          'circle-tuple Bezout criterion; fixed-parameter exclusion identities; '
          +str(cases)+' integer target checks.')


if __name__ == '__main__':
    main()
