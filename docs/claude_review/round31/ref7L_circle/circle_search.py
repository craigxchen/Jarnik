"""Search for five marks with ALL roots of N on the unit circle (r=0) and contact h, degree D.
Marks: mark j has root exp(i*sigma_jk*theta_k), k=1..D (theta_k in (0,2pi), with repetition allowed).
Contact >= h  <=>  p_m equal for m<h  <=>  sum_k (sigma_jk - sigma_1k) sin(m theta_k) = 0, m=1..h-1.
WLOG sigma_1 = (+...+) (global flips theta->-theta). Columns = 4-bit patterns of marks 2..5, column 0000 excluded
(common factor). (VII.4) with r=0 says 20h <= 9D.  Usage: python3 circle_search.py D h starts"""
import sys, itertools, numpy as np
from scipy.optimize import least_squares
D, h, starts = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
part = int(sys.argv[4]) if len(sys.argv) > 4 else 0; nparts = int(sys.argv[5]) if len(sys.argv) > 5 else 1
types = list(range(1, 16))
rng = np.random.default_rng(1 + part)
def sig(code):
    S = np.ones((5, D))
    for k, c in enumerate(code):
        for j in range(4):
            if (c >> j) & 1: S[j+1, k] = -1
    return S
def true_contact(S, th):
    # compute actual contact numerically from power sums of roots
    Z = np.exp(1j*S*th[None, :])
    for m in range(1, 3*D):
        p = (Z**m).sum(axis=1)
        if np.max(np.abs(p - p[0])) > 1e-7: return m
    return 99
best = []
count = 0
for idx, code in enumerate(itertools.combinations_with_replacement(types, D)):
    if idx % nparts != part: continue
    S = sig(code)
    if len({tuple(r) for r in S}) < 5: continue
    Dl = S[1:] - S[0]
    if np.linalg.matrix_rank(Dl) < 3: continue
    count += 1
    def F(th):
        return np.concatenate([Dl @ np.sin(m*th) for m in range(1, h)])
    for s in range(starts):
        th0 = rng.uniform(0.05, 2*np.pi-0.05, D)
        sol = least_squares(F, th0, xtol=1e-15, ftol=1e-15, gtol=1e-15)
        if np.max(np.abs(sol.fun)) > 1e-11: continue
        th = np.mod(sol.x, 2*np.pi)
        # nondegeneracy: no real roots, marks distinct as root multisets
        if np.min(np.abs(np.sin(th))) < 1e-6: continue
        Z = np.exp(1j*S*th[None, :])
        polys = [np.poly(z) for z in Z]
        ok = all(np.max(np.abs(polys[a]-polys[b])) > 1e-6 for a in range(5) for b in range(a+1, 5))
        if not ok: continue
        # primitivity: no root common to all five marks
        allroots = [np.round(z, 6) for z in Z]
        common = set(allroots[0]).intersection(*[set(x) for x in allroots[1:]])
        if common: continue
        hc = true_contact(S, th)
        best.append((hc, code, th.round(6).tolist()))
        if 20*hc > 9*D: print('VIOLATION', hc, code, th.tolist(), flush=True)
        break
print(f'D={D} h>={h} part {part}: codes tried {count}, solutions {len(best)}; max contact found',
      max([b[0] for b in best], default=None), flush=True)
for b in sorted(best, reverse=True)[:5]: print(b)
