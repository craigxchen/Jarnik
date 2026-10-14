"""Referee: Loc_m(1) test in R = S/(x_i^3) for Sylvester n=8 boxes (M8 = L2, M8' = L3), own code.
usage: r_loc1.py L m1 m2 ..."""
import itertools, numpy as np, sys
from math import factorial
from flint import fmpz_mat, nmod_mat
L = int(sys.argv[1]); ms = [int(x) for x in sys.argv[2:]]
d = 7
cuts = [[1 if bin(a & (x + 1)).count('1') % 2 else 0 for x in range(d)] for a in range(1, 8)]  # rows x=1..7, base 0
Pm = np.array(cuts)
pts = np.array(list(itertools.product((-1, 0, 1), repeat=d)))
A = pts[(np.abs(pts @ Pm.T) <= L).all(1)]
print('L', L, '|A| =', len(A))
def mons(e, le=False): return [a for a in itertools.product(range(3), repeat=d) if (sum(a) <= e if le else sum(a) == e)]
Hc = []
for e in range(0, 15):
    rows = [[int(np.prod([k[x] ** a[x] for x in range(d)])) for k in A] for a in mons(e, True)]
    Hc.append(fmpz_mat(rows).rank())
    if Hc[-1] == len(A): break
D = [Hc[0]] + [Hc[i] - Hc[i - 1] for i in range(1, len(Hc))]
print('Delta H (exact):', D)
dirs = []
for v in itertools.product(range(-2, 3), repeat=d):
    v = np.array(v)
    if not v.any() or np.gcd.reduce(np.abs(v[v != 0])) != 1 or v[np.nonzero(v)[0][0]] < 0: continue
    dirs.append(v)
dirs = np.array(dirs); vals = A @ dirs.T
tau = np.array([len(np.unique(vals[:, i])) - 1 for i in range(len(dirs))])
p = 1000000007
for m in ms:
    Rm = mons(m); idx = {a: i for i, a in enumerate(Rm)}
    rows = []
    for v, t in zip(dirs, tau):
        if t + 1 > m: continue
        f = {}
        for c in itertools.combinations_with_replacement(range(d), t + 1):
            a = [0] * d
            for x in c: a[x] += 1
            a = tuple(a)
            if a in f or max(a) > 2: continue
            cf = factorial(t + 1)
            for x in range(d): cf = cf // factorial(a[x]) * int(v[x]) ** a[x]
            if cf: f[a] = cf
        for c in itertools.combinations_with_replacement(range(d), m - t - 1):
            h = [0] * d
            for x in c: h[x] += 1
            r = {}
            for a, cf in f.items():
                b = tuple(a[i] + h[i] for i in range(d))
                if max(b) <= 2: r[idx[b]] = (r.get(idx[b], 0) + cf) % p
            if r: rows.append(r)
    Mx = nmod_mat(max(1, len(rows)), len(Rm), [0] * (max(1, len(rows)) * len(Rm)), p)
    for i, r in enumerate(rows):
        for jj, c in r.items(): Mx[i, jj] = c
    tgt = len(Rm) - (D[m] if m < len(D) else 0)
    print('m=%d: dim R_m=%d, target dim J(A)_m/(x^3) = %d, rank_p(width powers) = %d -> %s' % (m, len(Rm), tgt, Mx.rank(), 'Loc_m(1) HOLDS (certified)' if Mx.rank() == tgt else 'deficit (mod p, evidence)'))
