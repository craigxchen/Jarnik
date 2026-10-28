from fractions import Fraction as Fr
from math import gcd, log
from itertools import combinations
from fam import family, fder, other_root, disc_T, fval
from ell import quotient_curve
from lift import lift_T
from profile import profile

def small_primes_of(n, bound=10**6):
    n = abs(n); ps = set(); p = 2
    while p * p <= n and p < bound:
        while n % p == 0:
            ps.add(p); n //= p
        p += 1
    if n > 1 and n < bound * bound:
        ps.add(n); n = 1
    return ps, n

xs = [Fr(2), Fr(3), Fr(-5), Fr(1)]
F, C, g, pts = quotient_curve(xs)
k = g[3]
maps = F['maps']
consts = []
consts += list(xs) + [x - y for x, y in combinations(xs, 2)] + [x + y for x, y in combinations(xs, 2)]
consts += [F['e'][0], F['e'][2] ** 2 - 4 * F['e'][1] * F['e'][3], F['lam']['b'], F['lam']['c'], F['lam']['b'] - F['lam']['c']]
Db, Dc = disc_T(maps['b']), disc_T(maps['c'])
s = -Db[1] / Db[2]; t = -Dc[0] / Dc[1]
consts += [s, t, s - t] + [x * x - s for x in xs] + [x * x - t for x in xs] + [x * x - y * y for x, y in combinations(xs, 2)]
for x in xs:
    T = x * x
    ds = [fder(maps[q], x) for q in 'abc']
    consts += ds + [u - v for u, v in combinations(ds, 2)]
    bp = other_root(maps['b'], T, x); cp = other_root(maps['c'], T, x)
    oth = [-x, bp, cp]
    consts += [u - v for u, v in combinations(oth, 2)] + [u - y for u in oth for y in xs]
    consts += [Db[2], Dc[1]] + [maps[q][1][k_] for q in 'bc' for k_ in range(3)] + [maps[q][0][k_] for q in 'bc' for k_ in range(3)]
# discriminant of E': Y^2 = X^3 + A2 X^2 + A4 X  -> disc = 16 A4^2 (A2^2 - 4 A4)
consts.append(16 * C.A4 ** 2 * (C.A2 ** 2 - 4 * C.A4))
Sigma = set()
for c in consts:
    c = Fr(c)
    if c == 0:
        print('zero constant!'); continue
    for part in (c.numerator, c.denominator):
        ps, rest = small_primes_of(part)
        Sigma |= ps
        if rest > 1 and rest >= 10**6:
            Sigma.add(('big', rest))
print('Sigma (small part):', sorted(p for p in Sigma if isinstance(p, int)))
Q1, Q2 = pts[0], pts[1]
Dl = C.add(Q2, (Q1[0], -Q1[1]))
R = Q1; observed = set()
for N in range(1, 15):
    R = C.add(R, Dl)
    L = lift_T(F, R[0] / k)
    try:
        d, b, Dm = profile(xs, L[0])
    except ZeroDivisionError:
        continue
    for A, B in combinations(d, 2):
        gg = gcd(d[A], d[B])
        if gg > 1:
            ps, rest = small_primes_of(gg)
            observed |= ps
            if rest > 1: observed.add(('rest', rest))
print('primes shared between two blocks (N<=14):', sorted(p for p in observed if isinstance(p, int)), [p for p in observed if not isinstance(p, int)])
print('all shared primes in Sigma:', all((p in Sigma) for p in observed))
