"""Exact direct-source identities; no synthetic large endpoint family asserted."""

from collections import Counter
from itertools import combinations
from math import gcd, isqrt, prod

from check_gaussian_reflection_replacement import ggcd, gnorm
from check_isometric_virtual_anchor import factor
from check_mobius_conductor_transfer import fixtures
from check_mobius_reciprocal_stretch_grid import realize


def primitive(points):
    common = points[0]
    for point in points[1:]:
        common = ggcd(common, point)
    return gnorm(common) == 1


def squareclass(value):
    return prod(p for p, e in factor(value).items() if e % 2)


def points_on_circle(N):
    radius = isqrt(N)
    result = set()
    for x in range(-radius, radius + 1):
        y = isqrt(N - x*x)
        if x*x + y*y == N:
            result.update(((x, y), (x, -y)))
    return sorted(result)


def audit(points):
    N = gnorm(points[0])
    assert N % 2 == 1 and all(gnorm(p) == N for p in points)
    counts = [0, 0, 0, 0]
    for z0 in points:
        rows = []
        G = 2*N
        for z in points:
            if z == z0:
                continue
            dot = z0[0]*z[0] + z0[1]*z[1]
            cross = z0[0]*z[1] - z0[1]*z[0]
            if dot == -N:
                a, b = 0, 1
            else:
                common = gcd(N + dot, cross)
                a, b = (N + dot)//common, cross//common
            f = a*a + b*b
            assert gcd(a, b) == 1 and 2*N % f == 0
            d, c = 2*N//f, N - dot
            assert c > 0 and d*b*b == c and d*a*a == 2*N - c
            assert c == ((z[0]-z0[0])**2 + (z[1]-z0[1])**2)//2
            assert squareclass(c) == squareclass(d)
            G = gcd(G, c)
            rows.append((a, b, d, c))
            counts[0] += 1
        g0 = gcd(*z0)
        assert G in (g0, 2*g0)
        counts[1] += 1
        counts[3] += g0 > 1
        selected = list({row[3]: row for row in rows}.values())
        for size in (2, 4):
            for subset in combinations(selected, size):
                label_product = prod(row[2] for row in subset)
                root = isqrt(label_product)
                if root*root != label_product:
                    continue
                value = prod(2*N-row[3] for row in subset)
                assert value == (root*prod(row[0] for row in subset))**2
                normalized = prod((2*N-row[3])//G for row in subset)
                assert isqrt(normalized)**2 == normalized
                counts[2] += 1
        if N >= 4:
            multiplicities = Counter(squareclass(d) for a, b, d, c in rows
                                     if c*c <= 4*N)
            assert max(multiplicities.values(), default=0) <= 2
    return counts


def main():
    total = [0, 0, 0, 0]
    tuples = [realize(rows)[0] for rows in fixtures()]
    for N in (5, 13, 25, 65, 85, 125, 325, 625):
        points = points_on_circle(N)
        tuples.extend(combinations(points, 3))
    for points in tuples:
        if primitive(points):
            result = audit(points)
            total = [a+b for a, b in zip(total, result)]
    assert total[3] > 0
    print(f"PASS: {total[0]} primitive pair dictionaries; {total[1]} exact "
          f"anchor gcds ({total[3]} nonprimitive anchors); {total[2]} "
          "even-product squares; all radial-cap squareclass counts.")


if __name__ == '__main__':
    main()
