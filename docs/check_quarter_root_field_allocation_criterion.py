"""Exact source-exponent and actual Gaussian quarter-root field checks."""

from itertools import combinations, product


POSITIVE = list(combinations(range(1, 8), 3))


def characters(a):
    total = sum(a)
    return [2 * (a[0] + sum(a[i] for i in I)) - total for I in POSITIVE]


def mul(z, w):
    return z[0] * w[0] - z[1] * w[1], z[0] * w[1] + z[1] * w[0]


def power(z, n):
    value = (1, 0)
    for _ in range(n):
        value = mul(value, z)
    return value


def norm(z):
    return z[0] ** 2 + z[1] ** 2


if __name__ == "__main__":
    count = 0
    orientations = []
    for a in product(range(3), repeat=8):
        if sum(a) % 2 == 0:
            continue
        values = characters(a)
        b = sorted(a[1:])
        positive = a[0] > sum(b[3:]) - sum(b[:3])
        negative = a[0] < sum(b[:4]) - sum(b[4:])
        assert positive == all(v > 0 for v in values)
        assert negative == all(v < 0 for v in values)
        assert (len({v % 4 for v in values}) == 1) == (len({v % 2 for v in b}) == 1)
        if len(orientations) < 64:
            orientations.append(tuple(int(v < 0) for v in values))
        count += 1
    for f in orientations:
        for g in orientations:
            # Four valuation coordinates: pi_p, bar(pi_p), pi_q, bar(pi_q).
            classes = [(1 << u) | (1 << (2 + v)) for u, v in zip(f, g)]
            target = classes[0]
            can_conjugate_to_one_class = all(c in (target, target ^ 15) for c in classes)
            assert can_conjugate_to_one_class == (len({u ^ v for u, v in zip(f, g)}) == 1)
    pi, barpi = (-1, 2), (-1, -2)
    for e in range(16, 41):
        a = (e, 0, 1, 2, 3, 4, 5, e - 10)
        rows = [mul(power(pi, v), power(barpi, e - v)) for v in a]
        assert len(set(rows)) == 8
        assert all(norm(z) == 5**e for z in rows)
        assert min(a) == 0 and max(a) == e
        values = characters(a)
        assert min(values) == 1 and all(v % 2 == 1 for v in values)
        assert {v % 4 for v in values} == {1, 3}
        for value in values:
            A = power(pi, value)
            h = 5 ** ((value - 1) // 2)
            beta = power(pi, (value - 1) // 2)
            assert norm(A) == 5 * h**2
            assert A == mul(pi, mul(beta, beta))
        v = [0] * 5
        for ell in range(1, e + 1):
            r = sum(x >= ell for x in a)
            v[min(r, 8 - r)] += 1
        assert v[1:] == [11, e - 14, 2, 1]
        assert 9 * v[1] + 4 * v[2] + v[3] - 2 * e == 2 * e + 45
    print(f"Passed {count} odd-total allocation criteria and 25 actual Gaussian families.")
    print("Passed 4096 two-prime independent-conjugation squareclass checks.")
