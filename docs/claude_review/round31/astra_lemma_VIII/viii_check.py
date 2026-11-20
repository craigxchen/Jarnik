"""Exact checks for the two integer matrix encodings in Part VIII."""

from itertools import combinations, product
from math import gcd, prod


def factor(n):
    p = 2
    while p * p <= n:
        e = 0
        while n % p == 0:
            n //= p
            e += 1
        if e:
            yield p, e
        p += 1
    if n > 1:
        yield n, 1


def verify(columns, blocks, cofactors):
    assert columns[0] == (1, 0)
    assert all(gcd(x, y) == 1 for x, y in columns)
    for i, j in combinations(range(len(columns)), 2):
        x, y = columns[i]
        u, v = columns[j]
        forced = prod(n for side, n in blocks if i in side and j in side)
        assert x * v - y * u == cofactors(i, j) * forced


def crt_sample():
    # Three odd primes, one 2-power, one constant prime and coalescing classes.
    moduli = [16, 27, 25, 7, 11]
    residues = [(1, 9), (5, 14), (9, 12), (2, 2), (6, 10)]
    q = prod(moduli)
    terms = [q // v * pow(q // v, -1, v) for v in moduli]
    return q, sorted({sum(a * c for a, c in zip(row, terms)) % q
                      for row in product(*residues)})


def congruence_data(q, g, xs):
    ys = [(x - xs[0]) // g for x in xs]
    Q = q // g
    ds = [gcd(y, Q) for y in ys]
    columns = [(Q // d, y // d) for y, d in zip(ys, ds)]
    blocks = []
    for p, e in factor(Q):
        pe = p**e
        classes = {y % pe for y in ys}
        assert len(classes) <= 2 and all(gcd(c, p) == 1 for c in classes if c)
        blocks.append(({i for i, y in enumerate(ys) if y % pe}, pe**2))

    def cofactor(i, j):
        difference = ys[j] - ys[i]
        return difference // gcd(Q, difference)

    verify(columns, blocks, cofactor)


def divisor_data(n, q, xs):
    a = xs[0]
    gs = [gcd(a, x) for x in xs]
    columns = [(a // g, (x - a) // (q * g)) for x, g in zip(xs, gs)]
    blocks = []
    for p, e in factor(n):
        for s in range(1, e + 1):
            pe = p**s
            side = {i for i, x in enumerate(xs) if (x % pe == 0) != (a % pe == 0)}
            blocks.append((side, p))

    def cofactor(i, j):
        v = q * gcd(xs[i], xs[j])
        assert (xs[j] - xs[i]) % v == 0
        return (xs[j] - xs[i]) // v

    verify(columns, blocks, cofactor)


def square_reduction(q, a, roots):
    g = Q = 1
    for p, e in factor(q):
        v = 0
        while v < e and a % p**(v + 1) == 0:
            v += 1
        if v == e:
            g *= p**((e + 1) // 2)
        else:
            assert v % 2 == 0
            g *= p**(v // 2)
            Q *= p**(e - v)
    assert q <= g * g * Q and all(x % g == 0 for x in roots)
    Q2 = Q if Q % 2 else Q // 2
    for p, e in factor(Q2):
        assert len({(x // g) % p**e for x in roots}) <= 2


def main():
    q, marks = crt_sample()
    common_divisor = gcd(q, *(x - marks[0] for x in marks))
    n_crt = 0
    for eight in combinations(marks, 8):
        congruence_data(q, common_divisor, eight)
        n_crt += 1

    n = 2**4 * 3**3 * 5**2 * 7**2
    divisors = sorted(prod(t) for t in product(*[[p**s for s in range(e + 1)]
                                                  for p, e in factor(n)]))
    n_div = 0
    for q in range(1, 17):
        for r in range(q):
            if gcd(r, q) != 1:
                continue
            xs = [d for d in divisors if d % q == r]
            if len(xs) >= 8:
                for start in range(len(xs)):
                    selected = (xs[start:] + xs[:start])[:8]
                    divisor_data(n, q, selected)
                    n_div += 1

    n_sq = 0
    for q in range(1, 301):
        by_residue = {}
        for x in range(q):
            by_residue.setdefault(x * x % q, []).append(x)
        for a, roots in by_residue.items():
            square_reduction(q, a, roots)
            n_sq += 1
    print(f"8-column CRT matrices: {n_crt}; divisor matrices: {n_div}; "
          f"nonempty square-root classes: {n_sq}")


if __name__ == "__main__":
    main()
