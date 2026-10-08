"""Exact arithmetic audits for the one-sided cone and angular-collapse note."""

from functools import reduce
from itertools import combinations
from math import atan, gcd, lcm, log, prod, sqrt


def mul(z, w):
    return z[0] * w[0] - z[1] * w[1], z[0] * w[1] + z[1] * w[0]


def conj(z):
    return z[0], -z[1]


def norm(z):
    return z[0] ** 2 + z[1] ** 2


def power(z, n):
    out = (1, 0)
    for _ in range(n):
        out = mul(out, z)
    return out


def ggcd(z, w):
    while w != (0, 0):
        numerator = mul(z, conj(w))
        denominator = norm(w)
        q = tuple((2 * n + denominator) // (2 * denominator)
                  for n in numerator)
        product = mul(q, w)
        z, w = w, (z[0] - product[0], z[1] - product[1])
    return z


def realize(s, d, multiplier=1):
    r = len(d)
    assert len({tuple(row) for row in s}) == len(s)
    assert all({row[j] for row in s} == {-1, 1} for j in range(r))
    a = [2 * (j + 1) for j in range(r)]
    scale = lcm(*(abs(a[j] ** 2 - a[k] ** 2)
                  for j, k in combinations(range(r), 2))) if r > 1 else 1
    t = scale * multiplier
    blocks = [(x * t, 1) for x in a]
    assert all(norm(ggcd(z, conj(z))) == 1 for z in blocks)
    assert all(gcd(norm(z), norm(w)) == 1
               for z, w in combinations(blocks, 2))
    gammas = [power(z, e) for z, e in zip(blocks, d)]
    points = [reduce(mul, (z if sign == 1 else conj(z)
                           for z, sign in zip(gammas, row)), (1, 0))
              for row in s]
    n = prod(norm(z) for z in gammas)
    assert all(norm(z) == n for z in points)
    assert len(set(points)) == len(points)
    assert norm(reduce(ggcd, points)) == 1
    phi = [e * atan(1 / (x * t)) for e, x in zip(d, a)]
    theta = [sum(sign * p for sign, p in zip(row, phi)) for row in s]
    delta = max(theta) - min(theta)
    assert delta <= 2 * sum(e / (x * t) for e, x in zip(d, a))
    if s[0] == [1] * r:
        assert sum(phi) < 1.5707963267948966
        largest = 2 * max(sum(p for sign, p in zip(row, phi) if sign == -1)
                          for row in s)
        assert abs(delta - largest) <= 1e-12 * largest
        # These checks are numerical diagnostics of the proved inequalities.
        assert log(delta / 2) >= -log(n) / (2 * r) - 1e-12
        for p, z in zip(phi, gammas):
            assert z[1] != 0
            assert abs(z[1]) >= 1
            assert log(norm(z)) / 2 + log(delta / 2) >= -1e-12
    return blocks, points, n, t, delta


def paley(q):
    squares = {x * x % q for x in range(1, q)}

    def chi(x):
        x %= q
        return 0 if x == 0 else (1 if x in squares else -1)

    return [[1] * q] + [[-1 if x == j else -chi(x - j)
                         for j in range(q)] for x in range(q)]


def check_paley(q):
    s = paley(q)
    assert all(sum(a * b for a, b in zip(x, y)) == -1
               for x, y in combinations(s, 2))
    masks = [sum((entry == -1) << j for j, entry in enumerate(row))
             for row in s]
    count = 0
    largest_t = -q
    for rows in combinations(range(q + 1), 4):
        mask = 0
        for x in rows:
            mask ^= masks[x]
        tq = q - 2 * bin(mask).count("1")
        # Exact version of |tq| <= 3 sqrt(q)+4.
        excess = max(0, abs(tq) - 4)
        assert excess * excess <= 9 * q
        largest_t = max(largest_t, tq)
        count += 1
    blocks, points, n, t, delta = realize(s, [1] * q)
    weights = [log(norm(z)) for z in blocks]
    assert all(sum(w * a * b for w, a, b in zip(weights, x, y)) < -log(4)
               for x, y in combinations(s, 2))
    # The all-plus row is the maximum angle, with no possible phase cancellation.
    assert points[0] == reduce(mul, blocks, (1, 0))
    assert log(delta) + log(n) / 4 > 0
    print(f"q={q}: {count} exact quadruples; min limiting defect "
          f"{q-largest_t}/{q}; actual primitive angular-collapse tuple checked")


def check_weighted_general_profiles():
    fixtures = [
        ([[1, 1, 1], [-1, 1, 1], [1, -1, 1], [1, 1, -1]], [1, 2, 3]),
        ([[1, -1, 1, -1], [-1, 1, -1, 1],
          [1, 1, -1, -1], [-1, -1, 1, 1]], [3, 1, 2, 4]),
        ([[1, 1], [-1, 1], [1, -1], [-1, -1]], [1, 1]),
        ([[1], [-1]], [2]),
    ]
    for s, d in fixtures:
        for multiplier in (1, 10, 100):
            realize(s, d, multiplier)
    for r in range(3, 100):
        assert 2 * r / (r - 2) <= 6
    print("12 weighted actual fixtures; r=1,2 boundaries; uniform exponent checked")


if __name__ == "__main__":
    check_weighted_general_profiles()
    for q in (11, 19, 31, 43):
        check_paley(q)
    print("PASS: exact Gaussian realizations and profile checks; "
          "phase floating-point values are diagnostics only")
