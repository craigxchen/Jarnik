"""Exact certificate, nested-prime, and pairing checks for isolable patterns."""

from fractions import Fraction
from itertools import combinations
from math import isqrt


def mul(z, w):
    a, b = z
    c, d = w
    return a * c - b * d, a * d + b * c


def conj(z):
    return z[0], -z[1]


def power(z, n):
    out = (1, 0)
    for _ in range(n):
        out = mul(out, z)
    return out


def split_primes(n):
    out = []
    p = 5
    while len(out) < n:
        if p % 4 == 1 and all(p % d for d in range(2, isqrt(p) + 1)):
            out.append(p)
        p += 2
    return out


def gaussian_prime(p):
    for a in range(1, isqrt(p) + 1):
        b2 = p - a * a
        b = isqrt(b2)
        if b and b * b == b2:
            return a, b
    raise AssertionError(p)


def invert(a):
    n = len(a)
    m = [[Fraction(a[i][j]) for j in range(n)]
         + [Fraction(i == j) for j in range(n)] for i in range(n)]
    for j in range(n):
        pivot = next((i for i in range(j, n) if m[i][j]), None)
        if pivot is None:
            return None
        m[j], m[pivot] = m[pivot], m[j]
        scale = m[j][j]
        m[j] = [x / scale for x in m[j]]
        for i in range(n):
            if i != j:
                scale = m[i][j]
                m[i] = [x - scale * y for x, y in zip(m[i], m[j])]
    return [row[n:] for row in m]


def sign(cut, i):
    return 1 if i in cut else -1


def check_nested_exact_identity():
    m, k = 10, 5
    a = list(range(k))
    b = list(range(k, m))
    regular = [set((i, j)) for i in a for j in b]
    nested_big = {0, 1, 5, 6}
    nested_small = {0, 5}
    star_cut = set(a)
    cuts = regular + [nested_big, nested_small, star_cut]
    ps = split_primes(len(regular) + 2)
    pis = [gaussian_prime(p) for p in ps]
    common = (2, 1)
    units = [(1, 0), (0, 1), (-1, 0), (0, -1)]
    star = power(pis[-1], 3)
    rows = []
    for i in range(m):
        z = mul(units[i % 4], common)
        for cut, pi in zip(regular, pis[:len(regular)]):
            z = mul(z, pi if i in cut else conj(pi))
        # One rational prime has two genuinely nested regular layers.
        ap = int(i in nested_big) + int(i in nested_small)
        z = mul(z, mul(power(pis[-2], ap), power(conj(pis[-2]), 2 - ap)))
        z = mul(z, star if i in star_cut else conj(star))
        rows.append(z)
    lam = [1 if i in a else -1 for i in range(m)]
    assert sum(lam) == 0
    for cut in cuts[:-1]:
        assert sum(lam[i] * sign(cut, i) for i in range(m)) == 0
    assert sum(lam[i] * sign(star_cut, i) for i in range(m)) == m
    prod_a = prod_b = u = (1, 0)
    for i in a:
        prod_a = mul(prod_a, rows[i])
        u = mul(u, units[i % 4])
    for j in b:
        prod_b = mul(prod_b, rows[j])
        u = mul(u, conj(units[j % 4]))
    assert mul(prod_a, power(conj(star), k)) == mul(mul(u, prod_b), power(star, k))
    assert len(set(rows)) == m
    return len(cuts)


def check_pairing():
    m = 12
    examples = []
    base = set(range(6))
    choices = [base, {0, 1, 2, 3, 4, 6}, {0, 1, 2, 3, 7, 8},
               {0, 1, 2, 5, 9, 10}]
    for q in (2, 3, 4):
        columns = choices[:q]
        e = [[sign(cut, i) for i in range(m)] for cut in columns]
        assert all(sum(col) == 0 for col in e)
        h = [[sum(x * y for x, y in zip(e[a], e[b])) for b in range(q)]
             for a in range(q)]
        inv = invert(h)
        assert inv is not None
        for a in range(q):
            for b in range(q):
                assert sum(inv[a][c] * h[c][b] for c in range(q)) == (a == b)
        costs = [Fraction(m, 2) * sum(abs(x) for x in inv[j]) for j in range(q)]
        assert all(cost >= Fraction(1, 2) for cost in costs)
        assert any(cost > Fraction(1, 2) for cost in costs)
        examples.append((q, max(costs)))

    # Four distinct balanced columns can have a genuine rank-three relation.
    e1 = (1, 1, 1, 1, -1, -1, -1, -1)
    e2 = (1, 1, -1, -1, 1, 1, -1, -1)
    e3 = (1, 1, 1, -1, -1, 1, -1, -1)
    e4 = tuple(a + b - c for a, b, c in zip(e1, e2, e3))
    es = (e1, e2, e3, e4)
    assert all(set(col) == {-1, 1} and sum(col) == 0 for col in es)
    assert len(set(es)) == 4
    h = [[sum(x * y for x, y in zip(es[a], es[b])) for b in range(4)]
         for a in range(4)]
    assert invert(h) is None
    return examples


def check_centered_pairing():
    m = 9
    cuts = [{0, 1, 2, 3}, {0, 1, 4, 5, 6}]
    e = [[Fraction(sign(cut, i)) for i in range(m)] for cut in cuts]
    c = [[x - sum(col) / m for x in col] for col in e]
    assert any(sum(col) for col in e)
    assert all(sum(col) == 0 for col in c)
    k = [[sum(x * y for x, y in zip(c[a], e[b])) for b in range(2)]
         for a in range(2)]
    inv = invert(k)
    assert inv is not None
    for j in range(2):
        lam = [sum(inv[j][a] * c[a][i] for a in range(2)) for i in range(m)]
        assert sum(lam) == 0
        assert [sum(x * y for x, y in zip(lam, e[b])) for b in range(2)] == \
               [Fraction(j == b) for b in range(2)]
    costs = [sum(abs(inv[j][a]) * sum(abs(x) for x in c[a])
                 for a in range(2)) / 2 for j in range(2)]
    assert all(cost > 0 for cost in costs)
    return costs


if __name__ == "__main__":
    columns, costs = check_nested_exact_identity(), check_pairing()
    centered = check_centered_pairing()
    print(f"PASS: {columns} literal source columns with two nested regular layers, "
          f"exact arbitrary-unit star isolation; q=2,3,4 nonorthogonal inverses "
          f"and costs {costs}; centered unbalanced q=2 costs {centered}; "
          "a rank-deficient four-column witness.")
