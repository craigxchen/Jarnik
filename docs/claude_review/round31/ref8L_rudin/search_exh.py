"""Exhaustive search, small diameters: sets X of B integers with min 0 and max L,
Qmax(X) = prod_p p^{e_p}, e_p = largest e such that X meets at most two
classes mod p^e.  Exponent E(X) = log Qmax / log L.  VIII.2 asserts: for some
B0, eps, every X with |X| > B0 and Qmax large has E(X) < 2/(1+2 eps) < 2.
CG Theorem 9 (pigeonhole) gives E(X) <= binom(B,2)/m(B), m(B)=min same-class pairs.
Usage: search_exh.py B Lmax
"""
import sys
from math import log, comb
from itertools import combinations

B = int(sys.argv[1])
LMAX = int(sys.argv[2])

spf = list(range(LMAX + 2))
for i in range(2, int((LMAX + 1) ** 0.5) + 1):
    if spf[i] == i:
        for j in range(i * i, LMAX + 2, i):
            if spf[j] == j:
                spf[j] = i


def primes_of(n):
    s = set()
    while n > 1:
        p = spf[n]
        s.add(p)
        while n % p == 0:
            n //= p
    return s


def qmax(X):
    cand = primes_of(X[1] - X[0]) | primes_of(X[2] - X[0]) | primes_of(X[2] - X[1])
    Q = 1
    for p in cand:
        pe = p
        e = 0
        while len({x % pe for x in X}) <= 2:
            e += 1
            pe *= p
        Q *= p ** e
    return Q


def bound(B):
    a = B // 2
    m = comb(a, 2) + comb(B - a, 2)
    return comb(B, 2) / m


best = {}
overall = (0, None)
for L in range(B - 1, LMAX + 1):
    bL = (0, None)
    for mid in combinations(range(1, L), B - 2):
        X = (0,) + mid + (L,)
        Q = qmax(X)
        if Q > bL[0]:
            bL = (Q, X)
    E = log(bL[0]) / log(L) if L > 1 and bL[0] > 1 else 0
    best[L] = (E, bL[0], bL[1])
    if E > overall[0]:
        overall = (E, L, bL[0], bL[1])
print(f"B={B} Lmax={LMAX} pigeonhole bound={bound(B):.4f}")
print("overall best exponent:", overall)
for L in sorted(best):
    if L >= LMAX - 4 or L % max(1, LMAX // 8) == 0:
        print(L, "%.4f" % best[L][0], best[L][1], best[L][2])
