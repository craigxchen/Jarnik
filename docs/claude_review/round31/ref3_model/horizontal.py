"""Independent implementation of the III.6 horizontal CRT construction (columns (X_i,1))."""
import sys, random
from math import gcd, prod, isqrt
from itertools import combinations
from fractions import Fraction
from sympy import isprime, sqrt_mod, legendre_symbol
from sympy.ntheory.modular import crt
from mpmath import mp, mpf, atan2
m = int(sys.argv[1]); P0 = int(sys.argv[2])
assert isprime(m) and m % 4 == 3
chi = lambda x: 0 if x % m == 0 else legendre_symbol(x % m, m)
# Paley double core: columns (a,eps); entry row i: chi(a-i) for i!=a; at i=a: +1 for eps=+ (first family), -1 for eps=-
cols = [(a, e) for e in (1, -1) for a in range(m)]
def entry(i, c):
    a, e = c
    return e if i == a else chi(a - i)
A = [[1 if entry(i, c) == 1 else 0 for c in cols] for i in range(m)]
first = [k for k, c in enumerate(cols) if c[1] == 1]
for i in range(m):
    assert sum(A[i]) == m and sum(A[i][k] for k in first) == (m+1)//2
for i, j in combinations(range(m), 2):
    assert sum(A[i][k]*A[j][k] for k in range(2*m)) == (m-1)//2
    assert sum(A[i][k]*A[j][k] for k in first) == (m+1)//4
    assert sum(A[i][k] != A[j][k] for k in range(2*m)) == m+1
print("profile checks: row ones m, first-family ones (m+1)/2, shared (m-1)/2, first-family shared (m+1)/4, separation m+1")
# primes: 2m consecutive split primes above P0, random assignment
ps = []; q = P0
while len(ps) < 2*m:
    q += 1
    if q % 4 == 1 and isprime(q): ps.append(q)
rng = random.Random(2026); rng.shuffle(ps)
P = min(ps)
n = [prod(ps[k] for k in range(2*m) if A[i][k]) for i in range(m)]
g = {(i,j): prod(ps[k] for k in range(2*m) if A[i][k] and A[j][k]) for i,j in combinations(range(m),2)}
N = prod(ps); mp.dps = len(str(N))*2+50
# delete first-family columns so every row loses >= 3 ones (greedy, deterministic)
deleted = []
while any(sum(A[i][k] for k in deleted) < 3 for i in range(m)):
    best = max((k for k in first if k not in deleted), key=lambda k: sum(A[i][k] for i in range(m) if sum(A[i][kk] for kk in deleted) < 3))
    deleted.append(best)
S = [k for k in first if k not in deleted]
print("deleted", len(deleted), "first-family columns; |S| =", len(S))
r = {k: int(sqrt_mod(-1, ps[k])) for k in S}
H0 = isqrt(P**m); G = P**((m-1)//2)
b = [2*m*m*i + i*i for i in range(m)]
X = []
for i in range(m):
    ones = [k for k in S if A[i][k]]
    zeros = [k for k in S if not A[i][k]]
    mods = [2] + [ps[k] for k in ones]; vals = [0] + [r[k] for k in ones]
    x0, Mi = crt(mods, vals); x0, Mi = int(x0), int(Mi)
    assert Mi <= 2*G // P**2 * 2
    tgt = H0 + G*b[i]
    base = tgt - (tgt - x0) % Mi
    for t in range(2*m+2):
        x = base + t*Mi
        f = x*x + 1
        if all(f % (ps[k]**2) != 0 for k in ones) and all(f % ps[k] != 0 for k in zeros):
            break
    else:
        raise SystemExit("no candidate")
    X.append(x)
# checks
for i in range(m):
    f = X[i]**2 + 1
    for k in S:
        p = ps[k]
        if A[i][k]: assert f % p == 0 and f % (p*p) != 0 and (X[i] - r[k]) % p == 0
        else: assert f % p != 0
print("exact prescribed Gaussian valuations on S (consistent factor X = r_p mod p): OK; columns (X_i,1) primitive trivially")
gs_ok = True; maxerr = 0; logratio = []; fullfail = 0
for i, j in combinations(range(m), 2):
    gS = prod(ps[k] for k in S if A[i][k] and A[j][k])
    d = X[j] - X[i]
    gs_ok &= (d % gS == 0)
    fullfail += (d % g[i,j] != 0)
    err = abs(Fraction(d, g[i,j]) - (b[j]-b[i]))
    maxerr = max(maxerr, err)
    logratio.append(float(mp.log(mpf(gS))/mp.log(mpf(g[i,j]))))
vals = sorted(abs(b[j]-b[i]) for i,j in combinations(range(m),2))
print("g_{ij,S} | det on all pairs:", gs_ok, "; full g_ij divides on", m*(m-1)//2 - fullfail, "pairs of", m*(m-1)//2)
print("log g_S / log g range: %.3f .. %.3f" % (min(logratio), max(logratio)))
print("max |det/g_ij - (b_j-b_i)| = %.4g ; distinct |b_j-b_i|: %s" % (float(maxerr), len(set(vals)) == len(vals)))
print("max | |V|^2/n_i - 1 | = %.4g" % max(abs(float(Fraction(X[i]**2+1, n[i]) - 1)) for i in range(m)))
ang = [2*atan2(1, mpf(x)) for x in X]
print("doubled arg diam * N^(1/4) = %s (needs <= 1/2)" % mp.nstr((max(ang)-min(ang))*mpf(N)**0.25, 6))
