# exact global minima of sigma by exhaustive enumeration (int64 numpy, exact integers), formula of Lemma 1.1
# (formula itself checked against the pattern definition in c1_definition.py)
import numpy as np, itertools, math, sys
from plib import chi_table, patterns, N_def
def circ(q):
    ch = chi_table(q); return np.array([[ch[(a - t) % q] for a in range(q)] for t in range(q)], dtype=np.int64)
def sigma_batch(L, q, C):
    f = L[:, :q]; F = f.sum(1, keepdims=True); h = F + f; c = f @ C
    return np.maximum(np.abs(h), np.abs(c)).sum(1) - q
def is_pair_multiple(L):
    nz = (L != 0).sum(1); return nz == 2   # sum zero & two nonzero entries => l(e_x - e_y)
def run_box(q, B):
    C = circ(q); n = q + 1; vals = np.arange(-B, B+1)
    best = None; bestlam = None; cnt = 0; pairmin = None
    # enumerate first n-1 coords, last (inf) = -sum, keep |inf| <= B
    for chunk in itertools.product(vals, repeat=3):
        rest = np.array(list(itertools.product(vals, repeat=n-1-3)), dtype=np.int64)
        L = np.concatenate([np.tile(np.array(chunk, dtype=np.int64), (len(rest), 1)), rest], 1)
        inf = -L.sum(1); keep = np.abs(inf) <= B
        L = np.concatenate([L[keep], inf[keep][:, None]], 1)
        L = L[(L != 0).any(1)]
        s = sigma_batch(L, q, C); pm = is_pair_multiple(L)
        if pm.any():
            m = s[pm].min(); pairmin = m if pairmin is None else min(pairmin, m)
        s2 = s[~pm]
        if len(s2):
            m = s2.min()
            if best is None or m < best: best = m; bestlam = L[~pm][s2.argmin()]
    return best, bestlam, pairmin
for q, B in [(7, 2), (11, 1)]:
    best, lam, pm = run_box(q, B)
    # verify argmin by definition
    P = patterns(q); sd = N_def(list(lam), P) - q
    print(f"q={q}, all lambda with |lambda_x|<={B}: min sigma over non-pair-multiples = {best} (argmin {lam.tolist()}, "
          f"definition gives {sd}); min over pair multiples = {pm}; (q-1)/2 = {(q-1)//2}, q/3 = {q/3:.3f}", flush=True)
# q = 7: ternary F=0 with n >= 6 (the case where Prop 2.5(3)'s displayed (q-3)/2 = 2 < q/3)
q = 7; C = circ(q)
L = np.array([v + (0,) for v in itertools.product([-1, 0, 1], repeat=q) if sum(v) == 0 and sum(map(abs, v)) >= 6], dtype=np.int64)
print(f"q=7 ternary F=0 n>=6: count {len(L)}, min sigma {sigma_batch(L, q, C).min()} (Prop 2.5(1): n-1 >= 5)", flush=True)
# q = 11: ternary F = 0 by n, and |F| = 1 by n; Prop 2.5 / 2.4 bounds
q = 11; C = circ(q)
allv = np.array(list(itertools.product([-1, 0, 1], repeat=q)), dtype=np.int64)
for Fv in [0, 1]:
    sel = allv[allv.sum(1) == Fv]
    L = np.concatenate([sel, -np.full((len(sel), 1), Fv)], 1)
    nn = np.abs(L).sum(1); s = sigma_batch(L, q, C)
    for n in sorted(set(nn.tolist())):
        if n <= 2: continue
        m = s[nn == n].min()
        if Fv == 0:
            bound = max(n - 1, (q*(n-2) + 4*n - n*n)/(n+2) if n <= (q+1)//2 else -1)
        else:
            bound = (q-1)/2 if n == 4 else (q+1)/2
        print(f"q=11 F={Fv} n={n}: min sigma {m}, proved bound {bound:.3f}, ok={m >= bound - 1e-9}", flush=True)
