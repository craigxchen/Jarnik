"""Referee: independent sdc7 (L=2) box; exact H_A(8); Loc_8(1) test (width powers, directions in {-2..2}^7);
uncertainty budget data (tau_0, T_core)."""
import itertools, numpy as np, time, sys
from flint import fmpz_mat, nmod_mat
q = 7; n = 8; L = 2
def leg(a): a %= q; return 0 if a == 0 else (1 if pow(a, 3, q) == 1 else -1)
# skew Paley: index 0 = inf, 1..7 = F_7 elements 0..6
S = np.zeros((n, n), dtype=int)
for x in range(n):
    for y in range(n):
        if x == y: continue
        if x == 0: S[x, y] = 1
        elif y == 0: S[x, y] = -1
        else: S[x, y] = leg((y - 1) - (x - 1))
K = np.eye(n, dtype=int) + S
assert (K + K.T == 2 * np.eye(n)).all() and (K @ K.T == n * np.eye(n)).all()
Dg = np.diag(K[:, 0]); K = Dg @ K @ Dg
assert (K[:, 0] == 1).all() and (K + K.T == 2 * np.eye(n)).all()
cuts = []
for a in range(1, n):
    Sa = {x for x in range(n) if K[x, a] == -1}
    for C in (Sa, Sa ^ {0, a}):
        P = frozenset(C) if 0 not in C else frozenset(range(n)) - frozenset(C)
        cuts.append(P)
cuts = sorted(set(cuts), key=sorted)
print('distinct cuts', len(cuts))
d = n - 1
pts = np.array(list(itertools.product((-1, 0, 1), repeat=d)))   # k' = (k_1..k_7), base 0
Pm = np.array([[1 if (x + 1) in P else 0 for x in range(d)] for P in cuts])
A = pts[(np.abs(pts @ Pm.T) <= L).all(1)]
print('|A| =', len(A))
# exact Hilbert function (cumulative) via fmpz ranks, cube-reduced monomials
def mons_le(e): return [a for a in itertools.product(range(3), repeat=d) if sum(a) <= e]
def evm(mons):
    return [[int(np.prod([k[x] ** a[x] for x in range(d)])) for k in A] for a in mons]
Hc = []
for e in range(0, 10):
    r = fmpz_mat(evm(mons_le(e))).rank(); Hc.append(r)
    if r == len(A): break
print('cumulative H (exact):', Hc, ' Delta:', [Hc[0]] + [Hc[i] - Hc[i - 1] for i in range(1, len(Hc))])
# widths
dirs = []
for v in itertools.product(range(-2, 3), repeat=d):
    v = np.array(v)
    if not v.any() or np.gcd.reduce(np.abs(v[v != 0])) != 1 or v[np.nonzero(v)[0][0]] < 0: continue
    dirs.append(v)
dirs = np.array(dirs)
vals = A @ dirs.T
tau = np.array([len(np.unique(vals[:, i])) - 1 for i in range(len(dirs))])
from collections import Counter
print('tau histogram (<=7):', sorted(Counter(tau[tau <= 7]).items()))
tau0 = len(np.unique(-A.sum(1))) - 1
print('tau_0 (base coordinate) =', tau0)
print('core column taus:', [len(np.unique(A @ Pm[i])) - 1 for i in range(len(cuts))])
# Loc_8(1): span of l_v^{tau+1} * S_{7-tau} in S_8 (all monomials of degree 8), mod p
p = 1000000007
m = 8
mon8 = [a for a in itertools.product(range(m + 1), repeat=d) if sum(a) == m]
idx = {a: i for i, a in enumerate(mon8)}
from math import comb, factorial
def powform(v, e):
    # (v.x)^e as dict
    out = {}
    for c in itertools.combinations_with_replacement(range(d), e):
        a = [0] * d
        for x in c: a[x] += 1
        a = tuple(a)
        if a in out: continue
        coef = factorial(e)
        for x in range(d): coef //= factorial(a[x])
        val = coef
        for x in range(d): val *= int(v[x]) ** a[x]
        if val: out[a] = val
    return out
rows = []
for v, t in zip(dirs, tau):
    if t + 1 > m: continue
    f = powform(v, t + 1)
    for c in itertools.combinations_with_replacement(range(d), m - t - 1):
        h = [0] * d
        for x in c: h[x] += 1
        r = {}
        for a, cf in f.items():
            b = tuple(a[i] + h[i] for i in range(d)); r[idx[b]] = (r.get(idx[b], 0) + cf) % p
        rows.append(r)
print('generator rows', len(rows), flush=True)
Mx = nmod_mat(len(rows), len(mon8), [0] * (len(rows) * len(mon8)), p)
for i, r in enumerate(rows):
    for jj, c in r.items(): Mx[i, jj] = c
rk = Mx.rank()
print('dim S_8 =', len(mon8), ' rank_p(width-power part of J_8) =', rk, ' target dim J(A)_8 =', len(mon8) - (Hc[8] - Hc[7]) if len(Hc) > 8 else '?')
