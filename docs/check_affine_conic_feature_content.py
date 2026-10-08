"""Exact checks of the affine conic identities; standard library only."""

from itertools import combinations
from math import gcd, isqrt


def content(values):
    result = 0
    for value in values:
        result = gcd(result, value)
    return result


def triangle(a, b, c):
    return (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0])


def cross(a, b):
    return (
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    )


def check(points):
    m = len(points)
    g = content(triangle(*(points[i] for i in indices)) for indices in combinations(range(m), 3))

    def p(i, j, k):
        return triangle(points[i], points[j], points[k]) // g

    D = p(0, 1, 2)
    rows = []
    uvw = []
    for i in range(3, m):
        u, v, w = p(0, i, 2), p(0, 1, i), p(1, 2, i)
        assert w == D - u - v
        uvw.append((u, v, w))
        rows.append((u * (u - D), 2 * u * v, v * (v - D)))

    d1 = tuple(points[1][j] - points[0][j] for j in range(2))
    d2 = tuple(points[2][j] - points[0][j] for j in range(2))
    gram = (sum(x*x for x in d1), sum(x*y for x, y in zip(d1, d2)), sum(x*x for x in d2))
    gram_content = content(gram)
    normal = tuple(x // gram_content for x in gram)
    discriminant = normal[0] * normal[2] - normal[1]**2
    assert discriminant > 0 and isqrt(discriminant)**2 == discriminant
    assert all(sum(x*y for x, y in zip(row, normal)) == 0 for row in rows)

    cross_products = []
    multipliers = []
    for i, j in combinations(range(m - 3), 2):
        u, v, w = uvw[i]
        U, V, W = uvw[j]
        A = v * V * p(2, i + 3, j + 3)
        C = u * U * p(1, i + 3, j + 3)
        E = w * W * p(0, i + 3, j + 3)
        cubic = (2 * A, A + C - E, 2 * C)
        value = cross(rows[i], rows[j])
        assert value == tuple(D * x for x in cubic)
        assert cubic[0] % normal[0] == 0
        multiplier = cubic[0] // normal[0]
        assert cubic == tuple(multiplier * x for x in normal)
        cross_products.append(value)
        multipliers.append(multiplier)
    h = content(multipliers)
    assert h > 0
    assert content(x for value in cross_products for x in value) == abs(D) * h
    covolume_squared = sum((x // h)**2 for x in multipliers)
    for component in range(3):
        values = [value[component] for value in cross_products]
        common = content(values)
        if common:
            assert sum((value // common)**2 for value in values) == covolume_squared

    for a, b, c, d, e, f in combinations(range(m), 6):
        assert p(a, b, c) * p(a, d, e) * p(b, d, f) * p(c, e, f) == p(a, b, d) * p(a, c, e) * p(b, c, f) * p(d, e, f)


def main():
    count = 0
    for norm in (25, 65, 85, 125, 325, 425, 1105):
        bound = isqrt(norm)
        points = [(x, y) for x in range(-bound, bound + 1) for y in range(-bound, bound + 1) if x*x + y*y == norm]
        for size in (6, 8, len(points)):
            if size > len(points):
                continue
            for shift in range(min(5, len(points))):
                chosen = (points[shift:] + points[:shift])[:size]
                check(chosen)
                count += 1
    print(f"PASS: {count} exact configurations; conic normals, cubic minors, content, relation covolumes, and all six-point identities.")


if __name__ == "__main__":
    main()
