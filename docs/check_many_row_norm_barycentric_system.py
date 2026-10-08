"""Exact certificates for many_row_norm_barycentric_system.md; stdlib only."""
from fractions import Fraction as F
from itertools import combinations
from math import gcd, prod
import random


def factor(n):
    out = {}
    d = 2
    while d * d <= n:
        while n % d == 0:
            out[d] = out.get(d, 0) + 1
            n //= d
        d += 1
    if n > 1:
        out[n] = out.get(n, 0) + 1
    return out


def poly_mul(a, b):
    out = [0] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            out[i + j] += x * y
    return out


def moments():
    rng = random.Random(98109)
    for m in range(4, 9):
        for _ in range(24):
            # H_full=2+i, r=2, g=5, s=1; every other block is a unit.
            nodes = []
            while len(nodes) < m:
                a, y = rng.randrange(-40, 41), rng.randrange(1, 16)
                x = 2 * y + 5 * a
                if gcd(x, y) != 1 or (x - y) % 2 == 0:
                    continue
                if any(a * yy == aa * y for aa, yy in nodes):
                    continue
                nodes.append((a, y))
            a, y = map(list, zip(*nodes))
            x = [2 * yy + 5 * aa for aa, yy in nodes]
            b = [F(aa, yy) for aa, yy in nodes]
            k = [(yy + 2 * aa) ** 2 + aa ** 2 for aa, yy in nodes]
            tau = [prod(a[i] * y[j] - a[j] * y[i]
                        for j in range(m) if j != i) for i in range(m)]
            c = [F(k[i], tau[i]) for i in range(m)]
            for j in range(m - 3):
                assert sum(c[i] * x[i] ** j * y[i] ** (m - 4 - j)
                           for i in range(m)) == 0
            mm = [sum(c[i] * y[i] ** (m - 4) * b[i] ** j
                      for i in range(m)) for j in range(m)]
            assert all(z == 0 for z in mm[:m - 3])
            e1 = sum(b)
            e2 = sum(b[i] * b[j] for i, j in combinations(range(m), 2))
            alpha = prod(y)
            u, v = mm[-3], mm[-2] - e1 * mm[-3]
            w = mm[-1] - e1 * mm[-2] + e2 * mm[-3]
            assert (alpha * u, alpha * v, alpha * w) == (5, 4, 1)
            assert alpha ** 2 * (4 * u * w - v ** 2) == 4
            for i in range(m):
                lhs = F(5 * a[i] ** 2 + 4 * a[i] * y[i] + y[i] ** 2,
                        y[i] ** 2)
                lhs /= prod(b[i] - b[j] for j in range(m) if j != i)
                assert lhs == alpha * c[i] * y[i] ** (m - 4)


def countermodel():
    ab = [(-5, 2), (-5, 3), (-2, 5), (3, 3),
          (-5, 1), (-4, 4), (-1, 5), (1, 5)]
    q = [(a * a + b * b, 2 * a, 1) for a, b in ab]
    pp = [poly_mul(q[i], q[i + 1]) for i in range(0, 8, 2)]
    assert all(pp[0][j] + pp[1][j] == pp[2][j] + pp[3][j]
               for j in range(5))
    ss = {2}
    for a, b in ab:
        ss.update(factor(b))
        ss.update(factor(a * a + b * b))
    for (c, b, _), (e, d, _) in combinations(q, 2):
        u, v = d - b, e - c
        resultant = c * u * u - b * u * v + v * v
        assert resultant != 0
        ss.update(factor(abs(resultant)))
    assert ss == {2, 3, 5, 7, 13, 17, 29, 37, 53, 73, 89, 101, 109}
    modulus = prod(p ** (1 + max(factor(c).get(p, 0) for c, _, _ in q))
                   for p in ss)
    assert modulus == 348442128341958211928640
    for h in range(1, 9):
        t = modulus * h
        norms = []
        for a, b in ab:
            c = a * a + b * b
            assert t % c == 0
            real, imag = 1 + t * a // c, -t * b // c
            nn = real * real + imag * imag
            assert nn == ((t + a) ** 2 + b * b) // c
            assert gcd(real, imag) == 1 and (real - imag) % 2 != 0
            assert all(nn % p != 0 for p in ss)
            norms.append(nn)
        assert all(gcd(u, v) == 1 for u, v in combinations(norms, 2))
        assert (986 * norms[0] * norms[1] + 522 * norms[2] * norms[3]
                == 832 * norms[4] * norms[5] + 676 * norms[6] * norms[7])
    edges = {(0, 1): -39, (0, 2): -1352, (0, 3): -12,
             (1, 2): -1, (1, 3): -87, (2, 3): -4}
    def edge(i, j):
        return edges[i, j] if i < j else -edges[j, i]
    tau = [prod(edge(i, j) for j in range(4) if j != i) for i in range(4)]
    kk = [221, 1, 1, 1]
    assert [F(kk[i], tau[i]) for i in range(4)] == [
        F(c, 2822976) for c in [-986, 832, -522, 676]]
    assert 10 * 10 + 11 * 11 == 221 and gcd(10, 11) == 1


def boolean_exponents():
    for m in range(4, 10):
        all_sets = [frozenset(c) for r in range(1, m + 1)
                    for c in combinations(range(m), r)]
        for i in range(m):
            total = 0
            for tt in all_sets:
                expected = int(tt == {i}) + (len(tt) - 2 if i not in tt and len(tt) >= 3 else 0)
                got = (2 - len(tt) if i in tt else 0) + max(0, len(tt) - 2)
                assert expected == got
                total += expected
            assert total == (m - 5) * 2 ** (m - 2) + m + 2
    m = 5
    for i in range(m):
        content = 0
        total_degree = 0
        for mask in range(1, 1 << m):
            tt = {j for j in range(m) if mask >> j & 1}
            r = len(tt)
            aexp = int(tt == {i}) + (r - 2 if i not in tt and r >= 3 else 0)
            h = aexp + int(i in tt) - int(r >= 3)
            hb = aexp
            expected = ((2, 1) if tt == {i} else
                        (1, 2) if r == 4 and i not in tt else
                        (1, 0) if r == 2 and i in tt else
                        (0, 1) if r == 3 and i not in tt else (0, 0))
            assert (h, hb) == expected
            content += min(h, hb)
            total_degree += h + hb
        assert content == 2 and total_degree == 14


if __name__ == '__main__':
    moments()
    countermodel()
    boolean_exponents()
    print('Exact norm moments, progression, coupled residuals, and five-row content verified.')
