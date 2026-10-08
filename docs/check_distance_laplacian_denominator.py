"""Exact checks of positive tree sums and circle distance minors."""

from fractions import Fraction as F
from itertools import combinations, product
from math import atan, gcd, lcm, prod


def determinant(matrix):
    a = [[F(v) for v in row] for row in matrix]
    out = F(1)
    for i in range(len(a)):
        pivot = next((j for j in range(i, len(a)) if a[j][i]), None)
        if pivot is None:
            return F(0)
        if pivot != i:
            a[i], a[pivot] = a[pivot], a[i]
            out = -out
        value = a[i][i]
        out *= value
        for j in range(i + 1, len(a)):
            ratio = a[j][i] / value
            for k in range(i + 1, len(a)):
                a[j][k] -= ratio * a[i][k]
    return out


def trees(n):
    for code in product(range(n), repeat=n - 2):
        degrees = [1 + code.count(i) for i in range(n)]
        edges = []
        for vertex in code:
            leaf = next(i for i in range(n) if degrees[i] == 1)
            edges.append((leaf, vertex))
            degrees[leaf] -= 1
            degrees[vertex] -= 1
        edges.append(tuple(i for i in range(n) if degrees[i] == 1))
        yield edges


def gaussian_mul(z, w):
    return z[0] * w[0] - z[1] * w[1], z[0] * w[1] + z[1] * w[0]


def triangle(a, b, c):
    return (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0])


def check(values):
    n = len(values)
    q = [a * a + 1 for a in values]
    assert all(a % 2 == 0 for a in values)
    assert all(gcd(q[i], q[j]) == 1 for i, j in combinations(range(n), 2))
    x = [(F(a * a - 1, a * a + 1), F(2 * a, a * a + 1)) for a in values]
    d = [[sum((x[i][c] - x[j][c]) ** 2 for c in (0, 1)) for j in range(n)] for i in range(n)]
    radius_squared = lcm(*(d[i][j].denominator for i, j in combinations(range(n), 2)))
    assert radius_squared == prod(q)
    lap = [[sum(d[i]) if i == j else -d[i][j] for j in range(n)] for i in range(n)]
    tau = sum(prod(d[i][j] for i, j in edges) for edges in trees(n))
    assert tau == determinant([row[:-1] for row in lap[:-1]])
    assert tau.denominator == radius_squared ** (n - 1)
    assert gcd(tau.numerator, radius_squared) == 1

    center = [sum(row[c] for row in x) / n for c in (0, 1)]
    b = [[sum((x[i][c] - center[c]) * (x[j][c] - center[c]) for c in (0, 1)) for j in range(n)] for i in range(n)]
    e2 = sum(b[i][i] * b[j][j] - b[i][j] ** 2 for i, j in combinations(range(n), 2))
    area_sum = F(0)
    for i, j, k in combinations(range(n), 3):
        area = triangle(x[i], x[j], x[k])
        area_sum += area ** 2
        assert area ** 2 == d[i][j] * d[i][k] * d[j][k] / 4
        assert determinant([[d[u][v] for v in (i, j, k)] for u in (i, j, k)]) == 8 * area ** 2
    assert n * e2 == area_sum
    for inds in combinations(range(n), 4):
        assert determinant([[d[i][j] for j in inds] for i in inds]) == 0
    c = list(range(n - 1)) + [-sum(range(n - 1))]
    lhs = sum(c[i] * c[j] * d[i][j] for i in range(n) for j in range(n))
    rhs = -2 * sum(sum(c[i] * x[i][coord] for i in range(n)) ** 2 for coord in (0, 1))
    assert lhs == rhs <= 0

    integral = []
    for i in range(n):
        z = (1, 0)
        for j, a in enumerate(values):
            z = gaussian_mul(z, (a, 1 if i == j else -1))
        integral.append(z)
        assert z[0] ** 2 + z[1] ** 2 == radius_squared
    integer_area_sum = sum(triangle(integral[i], integral[j], integral[k]) ** 2 for i, j, k in combinations(range(n), 3))
    assert n * radius_squared ** 2 * e2 == integer_area_sum
    return radius_squared


def main():
    examples = [[2, 4, 6], [2, 4, 6, 10], [2, 4, 6, 10, 14], [2, 4, 6, 10, 14, 16]]
    for t in (5, 10, 25, 100, 1000):
        values = [2 * (t + i) for i in range(4)]
        n = check(values)
        endpoint = 2 * (atan(1 / (2 * t)) - atan(1 / (2 * (t + 3)))) * n ** 0.25
        print(f"T={t}: four-point endpoint normalization {endpoint:.10f}; denominator(tau)=N^3 verified")
    for values in examples:
        check(values)
    print(f"Passed {len(examples) + 5} exact circle configurations, including six-point full-denominator checks.")


if __name__ == "__main__":
    main()
