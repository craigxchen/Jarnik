"""Exact (interval) computation of the one-prime gap set G_1(p) in [0, pi/2).

G_1(p) = { t : for all n >= 1, dist(n t, (pi/2)Z) >= arcsin(p^{-n/2}),
                          dist(n t, pi/4 + (pi/2)Z) >= arcsin(p^{-n/2}/sqrt 2) }.
Truncated at n <= nmax; the surviving set is a superset of G_1(p).
Uses mpmath with 60 digits.  Evidence (finite computation), plus an exact
check of which surviving components contain the genuine angles arg(a+bi), a^2+b^2=p.
"""
import sys
from mpmath import mp, mpf, asin, pi, sqrt, atan2
mp.dps = 60

def subtract(intervals, a, b):
    out = []
    for (x, y) in intervals:
        if b <= x or a >= y:
            out.append((x, y)); continue
        if x < a: out.append((x, a))
        if b < y: out.append((b, y))
    return out

def gapset(p, nmax):
    q = pi / 2
    I = [(mpf(0), q)]
    for n in range(1, nmax + 1):
        ea = asin(mpf(p) ** (-mpf(n) / 2))
        ed = asin(mpf(p) ** (-mpf(n) / 2) / sqrt(2))
        # n t in (k q - ea, k q + ea) or (k q + q/2 - ed, k q + q/2 + ed), t in [0,q)
        for k in range(0, n + 1):
            c = k * q
            I = subtract(I, (c - ea) / n, (c + ea) / n)
            c2 = k * q + q / 2
            I = subtract(I, (c2 - ed) / n, (c2 + ed) / n)
        if not I:
            return I, n
    return I, nmax

if __name__ == "__main__":
    p = int(sys.argv[1]); nmax = int(sys.argv[2])
    I, nstop = gapset(p, nmax)
    tot = sum(y - x for x, y in I)
    gens = []
    for a in range(1, p):
        for b in range(1, p):
            if a * a + b * b == p:
                gens.append(atan2(b, a))
    print(f"p={p} nmax={nmax} stopped_at={nstop} components={len(I)} measure={mp.nstr(tot, 8)} fraction={mp.nstr(tot/(pi/2), 8)}")
    for (x, y) in I[:20]:
        tags = [mp.nstr(g, 12) for g in gens if x - mpf(10)**-40 <= g <= y + mpf(10)**-40]
        print("  [", mp.nstr(x, 15), ",", mp.nstr(y, 15), "] len", mp.nstr(y - x, 5), "genuine:", tags)
