"""Judge check of X.1 (own code): for every squarefree q <= QMAX and every triple of positive
lengths a <= b <= c <= floor(CAP*sqrt q) that is locally admissible at every p | q
(all three 0 mod p, or exactly one 0 mod p and the other two with equal squares mod p),
check (1) (q g)^2 | H, (2) T = H/(qg)^2 is 0/3 mod 4 (Q=q/g odd) or 0 mod 4 (Q even),
(3) H != 0 implies max^2 >= q g.  Also record the smallest max^2/(q g) over nondegenerate
triples with g < q.  Exact integer arithmetic in numpy int64 (|H| <= 81 CAP^4 q^2 < 2^63)."""
import sys
from fractions import Fraction
from math import isqrt
import numpy as np

QMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 2000
CAP = float(sys.argv[2]) if len(sys.argv) > 2 else 3


def factor(n):
    f, p = [], 2
    while p * p <= n:
        if n % p == 0:
            e = 0
            while n % p == 0:
                n //= p; e += 1
            f.append((p, e))
        p += 1
    if n > 1:
        f.append((n, 1))
    return f


cache = {}


def triples(N):
    if N not in cache:
        r = np.arange(1, N + 1, dtype=np.int64)
        a, b, c = np.meshgrid(r, r, r, indexing="ij")
        m = (a <= b) & (b <= c)
        cache.clear()
        cache[N] = (a[m], b[m], c[m])
    return cache[N]


stats = dict(moduli=0, admissible=0, nondeg=0, fail_div=0, fail_mod4=0, fail_cap=0)
best = None
for q in range(2, QMAX + 1):
    f = factor(q)
    if any(e > 1 for _, e in f):
        continue
    stats["moduli"] += 1
    N = int(CAP * q ** 0.5)
    a, b, c = triples(N)
    ok = np.ones(a.shape, bool)
    g = np.ones(a.shape, np.int64)
    for p, _ in f:
        za, zb, zc = a % p == 0, b % p == 0, c % p == 0
        nz = za.astype(int) + zb + zc
        A, Bq, Cq = (a * a) % p, (b * b) % p, (c * c) % p
        one = (nz == 1) & ((za & (Bq == Cq)) | (zb & (A == Cq)) | (zc & (A == Bq)))
        three = nz == 3
        ok &= one | three
        g = np.where(three, g * p, g)
    a1, b1, c1, g1 = a[ok], b[ok], c[ok], g[ok]
    H = (a1 + b1 + c1) * (-a1 + b1 + c1) * (a1 - b1 + c1) * (a1 + b1 - c1)
    D = (q * g1) ** 2
    stats["admissible"] += len(H)
    bad = H % D != 0
    stats["fail_div"] += int(bad.sum())
    T = H // D
    Q = q // g1
    m4 = np.where(Q % 2 == 1, (T % 4 == 0) | (T % 4 == 3), T % 4 == 0)
    stats["fail_mod4"] += int((~m4).sum())
    nd = H != 0
    stats["nondeg"] += int(nd.sum())
    stats["fail_cap"] += int((nd & (c1 * c1 < q * g1)).sum())
    sel = nd & (g1 < q)
    for x, y, z, gg in zip(a1[sel], b1[sel], c1[sel], g1[sel]):
        r = Fraction(int(z) * int(z), q * int(gg))
        if best is None or r < best[0]:
            best = (r, q, int(gg), (int(x), int(y), int(z)))
print("QMAX", QMAX, "CAP", CAP, stats)
print("smallest max^2/(qg) over nondegenerate triples with g<q:", best)
