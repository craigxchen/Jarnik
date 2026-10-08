"""Exact Gaussian check of deletion factors and primitive Plücker columns."""

from itertools import combinations, product
from math import lcm
from random import Random

from check_quartet_matching_gcd_cut_budget import ggcd, matching_content, mul, norm, sub


def add(a, b):
    return a[0] + b[0], a[1] + b[1]


def divide(a, b):
    d = norm(b)
    numerator = mul(a, (b[0], -b[1]))
    assert numerator[0] % d == numerator[1] % d == 0
    return numerator[0] // d, numerator[1] // d


def xgcd(a, b):
    r0, r1 = a, b
    s0, s1, t0, t1 = (1, 0), (0, 0), (0, 0), (1, 0)
    while r1 != (0, 0):
        d = norm(r1)
        num = mul(r0, (r1[0], -r1[1]))
        q = (num[0] + d // 2) // d, (num[1] + d // 2) // d
        r0, r1 = r1, sub(r0, mul(q, r1))
        s0, s1 = s1, sub(s0, mul(q, s1))
        t0, t1 = t1, sub(t0, mul(q, t1))
    assert add(mul(s0, a), mul(t0, b)) == r0
    return r0, s0, t0


def gcd_list(values):
    result = (0, 0)
    for value in values:
        result = ggcd(result, value)
    return result


def sum_list(values):
    result = (0, 0)
    for value in values:
        result = add(result, value)
    return result


def audit(points):
    m = len(points)
    D = gcd_list(sub(z, points[0]) for z in points[1:])
    u = [divide(sub(z, points[0]), D) for z in points]
    assert norm(gcd_list(u)) == 1
    A_i = [gcd_list(sub(u[k], u[j]) for j, k in combinations(range(m), 2)
                    if j != i and k != i) for i in range(m)]
    assert all(norm(ggcd(a, b)) == 1 for a, b in combinations(A_i, 2))
    A = (1, 0)
    for a in A_i:
        A = mul(A, a)
    G = gcd_list(matching_content(q) for q in combinations(points, 4))
    assert norm(divide(G, mul(mul(D, D), A))) == 1
    edges = [[divide(mul(mul(A_i[i], A_i[j]), sub(u[j], u[i])), A)
              for j in range(m)] for i in range(m)]
    assert norm(gcd_list(edges[i][j] for i, j in combinations(range(m), 2))) == 1
    for i, j, k in combinations(range(m), 3):
        assert add(sub(mul(A_i[i], edges[j][k]), mul(A_i[j], edges[i][k])),
                   mul(A_i[k], edges[i][j])) == (0, 0)
    for i, j, k, l in combinations(range(m), 4):
        assert add(sub(mul(edges[i][j], edges[k][l]), mul(edges[i][k], edges[j][l])),
                   mul(edges[i][l], edges[j][k])) == (0, 0)
    coefficients, current = [], (0, 0)
    for a in A_i:
        current, left, right = xgcd(current, a)
        coefficients = [mul(left, c) for c in coefficients] + [right]
    assert norm(current) == 1
    coefficients = [mul(c, (current[0], -current[1])) for c in coefficients]
    assert sum_list(mul(c, a) for c, a in zip(coefficients, A_i)) == (1, 0)
    b_i = [sum_list(mul(coefficients[i], edges[i][j]) for i in range(m))
           for j in range(m)]
    assert all(norm(ggcd(a, b)) == 1 for a, b in zip(A_i, b_i))
    for i, j in combinations(range(m), 2):
        assert sub(mul(A_i[i], b_i[j]), mul(A_i[j], b_i[i])) == edges[i][j]
    U = sum_list(mul(mul(c, a), v) for c, a, v in zip(coefficients, A_i, u))
    for a, b, v in zip(A_i, b_i, u):
        assert mul(a, v) == add(mul(a, U), mul(A, b))
    return D, G, A_i


def main():
    rng = Random(20260916)
    pool = [(x, y) for x in range(-10, 11) for y in range(-10, 11)]
    for _ in range(600):
        audit(rng.sample(pool, rng.randrange(4, 10)))
    # Force genuine deletion content, rather than relying on random residues.
    for _ in range(100):
        factor = rng.randrange(2, 12), 1
        majority = rng.sample(pool, rng.randrange(3, 8))
        points = [mul(factor, z) for z in majority] + [(1, 0)]
        audit(points)
    for k, q in ((3, 1), (3, 2), (4, 1), (4, 2)):
        M = lcm(*(j * j - l * l for j in range(1, k + 1) for l in range(1, j)))
        x = [2 * j * M * q for j in range(1, k + 1)]
        points = []
        for signs in product((-1, 1), repeat=k):
            z = (1, 0)
            for real, sign in zip(x, signs):
                z = mul(z, (real, sign))
            points.append(z)
        assert len(set(points)) == 2 ** k
        assert len(set(map(norm, points))) == 1 and norm(gcd_list(points)) == 1
        D, G, A_i = audit(points)
        assert norm(D) == 4 and norm(G) == 16 and all(norm(a) == 1 for a in A_i)
    print("PASS: 600 arbitrary tuples, 100 forced-deletion tuples,")
    print("and four primitive 8/16-point circle cubes with exact Norm(G)=16.")


if __name__ == "__main__":
    main()
