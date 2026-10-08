"""Exact combinatorial and Gaussian checks for higher cut-kernel content.

This checks the coordinate-support minimum, the triangle/pair cut witnesses,
the binomial height identities, and two local split-prime fixtures. It does
not assert a spanning theorem for graph monomials in K_k.
"""

from itertools import combinations, product
from math import comb


def g(m: int, k: int, s: int) -> int:
    q = m // 2
    return max(0, 2 * s - m, k + s - q)


def monomial_support_checks() -> None:
    # Exhaust every multidegree-two monomial for small even m. A monomial
    # survives all |T|=q-k+1 collapse tests exactly when its P-support has
    # at least q+k labels.
    for m in (4, 6, 8, 10, 12):
        q = m // 2
        for k in range(0, m // 6 + 1):
            cut_size = q - k + 1
            for mask in range(1 << m):
                support = {i for i in range(m) if mask & (1 << i)}
                hit_all = all(support.intersection(T)
                              for T in combinations(range(m), cut_size))
                assert hit_all == (len(support) >= q + k)

            allowed = []
            for exponents in product(range(3), repeat=m):
                if sum(exponents) != m:
                    continue
                support = sum(a > 0 for a in exponents)
                if k == 0:
                    assert support >= q
                if support >= q + k:
                    allowed.append(exponents)

            assert allowed
            for s in range(m + 1):
                minimum = min(sum(a for a in row[:s]) for row in allowed)
                assert minimum == g(m, k, s), (m, k, s, minimum, g(m, k, s))

    print("Monomial support minima agree with g_k(s) for all tested rows and cuts.")


def graph_cut_checks() -> None:
    # Dynamic programming over 2k triangle components and q-3k doubled pairs.
    for m in range(4, 42, 2):
        q = m // 2
        for k in range(0, m // 6 + 1):
            costs = {0: 0}
            for _ in range(2 * k):
                updated = {}
                for size, value in costs.items():
                    for inside, extra in enumerate((0, 0, 1, 3)):
                        total = size + inside
                        updated[total] = min(updated.get(total, 10**9), value + extra)
                costs = updated
            for _ in range(q - 3 * k):
                updated = {}
                for size, value in costs.items():
                    for inside, extra in enumerate((0, 0, 2)):
                        total = size + inside
                        updated[total] = min(updated.get(total, 10**9), value + extra)
                costs = updated
            assert len(costs) == m + 1
            for s in range(m + 1):
                assert costs[s] == g(m, k, s), (m, k, s, costs[s], g(m, k, s))
            # The graph has q-k components, so an edge-free subset has at
            # most one vertex per component.
            assert 2 * k + (m - 6 * k) // 2 == q - k
    print("Triangle/pair graphs attain g_k(s) for every tested m,k,s.")


def height_identity_checks() -> None:
    for m in range(4, 42, 2):
        q = m // 2
        base_internal = sum(comb(m, s) * max(0, 2 * s - m) for s in range(m + 1))
        assert base_internal == q * comb(m, q)
        A0 = m * 2 ** (m - 2) - q * comb(m, q)
        for k in range(0, m // 6 + 1):
            delta = sum(comb(m, s) * max(0, k - abs(s - q)) for s in range(m + 1))
            weighted_g = sum(comb(m, s) * g(m, k, s) for s in range(m + 1))
            assert weighted_g == base_internal + delta
            Ak = m * 2 ** (m - 2) - weighted_g
            assert Ak == A0 - delta
    expected = {(6, 1): 16, (8, 1): 162, (10, 1): 1048,
                (12, 1): 5820, (12, 2): 3312, (18, 3): 357528}
    for (m, k), value in expected.items():
        A0 = m * 2 ** (m - 2) - (m // 2) * comb(m, m // 2)
        delta = sum(comb(m, s) * max(0, k - abs(s - m // 2)) for s in range(m + 1))
        assert A0 - delta == value
    print("Binomial gcd and primitive-height identities pass.")


def gmul(a, b):
    return (a[0] * b[0] - a[1] * b[1], a[0] * b[1] + a[1] * b[0])


def gconj(a):
    return (a[0], -a[1])


def gnorm(a):
    return a[0] * a[0] + a[1] * a[1]


def determinant(a, b):
    return a[0] * b[1] - a[1] * b[0]


def valuation(n: int, p: int) -> int:
    n = abs(n)
    assert n != 0, "the valuation of zero is not a finite integer"
    value = 0
    while n and n % p == 0:
        value += 1
        n //= p
    return value


def gaussian_rows(primes, orientation_sets, m):
    rows = []
    for i in range(m):
        row = (1, 0)
        for p, pi in primes.items():
            row = gmul(row, pi if i in orientation_sets[p] else gconj(pi))
        rows.append(row)
    return rows


def triangle_edges(triples):
    return [edge for a, b, c in triples for edge in ((a, b), (b, c), (c, a))]


def check_local_prime_fixtures() -> None:
    # m=6, k=1, cut {0,1,2,3}; two internal graph edges give g_1(4)=2.
    primes6 = {5: (1, 2), 13: (3, 2), 17: (1, 4)}
    sets6 = {13: {0, 1, 2, 3}, 5: {0, 2, 4}, 17: {0, 1, 4}}
    rows6 = gaussian_rows(primes6, sets6, 6)
    assert len(set(rows6)) == 6 and len({gnorm(z) for z in rows6}) == 1
    assert gnorm(rows6[0]) == 1105
    edges6 = triangle_edges(((0, 1, 4), (2, 3, 5)))
    vals6 = [valuation(determinant(rows6[i], rows6[j]), 13) for i, j in edges6]
    assert vals6 == [1, 0, 0, 1, 0, 0]
    assert sum(vals6) == g(6, 1, 4) == 2

    # m=12, k=2, cut {0,...,7}; four internal triangle edges give g_2(8)=4.
    primes12 = {13: (3, 2), 5: (1, 2), 17: (1, 4), 29: (2, 5)}
    sets12 = {
        13: set(range(8)),
        5: {0, 2, 4, 6, 8, 10},
        17: {0, 1, 4, 5, 8, 9},
        29: {0, 3, 5, 6, 9, 10},
    }
    rows12 = gaussian_rows(primes12, sets12, 12)
    assert len(set(rows12)) == 12 and len({gnorm(z) for z in rows12}) == 1
    assert gnorm(rows12[0]) == 32045
    edges12 = triangle_edges(((0, 1, 8), (2, 3, 9), (4, 5, 10), (6, 7, 11)))
    vals12 = [valuation(determinant(rows12[i], rows12[j]), 13) for i, j in edges12]
    assert vals12 == [1, 0, 0] * 4
    assert sum(vals12) == g(12, 2, 8) == 4
    print("Exact p=13 Gaussian fixtures attain the m=6,k=1 and m=12,k=2 cut exponents.")


if __name__ == "__main__":
    monomial_support_checks()
    graph_cut_checks()
    height_identity_checks()
    check_local_prime_fixtures()
    print("All higher-k evaluation-height checks passed.")
