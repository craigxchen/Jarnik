"""Independent LP upper bound on zeta(A) (largest sum of lengths of a lattice zonotope in A), T=1 boxes.
For functionals f (slab functionals k(S), coordinates k_x incl. the base coordinate), any zonotope
k0 + sum [0,L_i] v_i in A has sum_i L_i |f(v_i)| <= h_A(f).  If y >= 0 with sum_f y_f |f(v)| >= 1 for every
admissible generator v, then zeta <= sum_f y_f h_A(f).  Generators: all v in {-2..2}^d \ 0 (superset of the
nonzero differences of A), filtered by |f(v)| <= h_A(f).  Cutting-plane LP, final check over ALL v exactly
(rational y)."""
import sys, itertools, numpy as np
from fractions import Fraction as Fr
from scipy.optimize import linprog
from rlib import *
name, n, pats, L = sys.argv[1], int(sys.argv[2]), None, int(sys.argv[3])
pats = {'sdc11': lambda: sdc(11), 'pal11': lambda: pal(11), 'sdc7': lambda: sdc(7), 'pal7': lambda: pal(7)}[name]()
d = n - 1
A, others = box(n, pats, L, 0)
full = np.zeros((len(A), n), dtype=np.int64); full[:, others] = A; full[:, 0] = -A.sum(1)
funcs = [np.array([1 if x in S else 0 for x in range(n)]) for S in pats] + [np.eye(n, dtype=np.int64)[x] for x in range(n)]
Fm = np.array(funcs)                     # functionals on full coordinates
h = np.array([(full @ f).max() - (full @ f).min() for f in Fm])
print(name, 'L', L, '|A|', len(A), 'widths: slabs', sorted(set(h[:len(pats)].tolist())), 'coords', h[len(pats):].tolist(), flush=True)
def chunks():
    vals = np.array([-2, -1, 0, 1, 2], dtype=np.int64)
    head = np.array(list(itertools.product(vals, repeat=d-6)), dtype=np.int64)
    tail = np.array(list(itertools.product(vals, repeat=6)), dtype=np.int64)
    for hrow in head:
        V = np.hstack([np.repeat(hrow[None, :], len(tail), 0), tail])
        Vf = np.hstack([-V.sum(1, keepdims=True), V])      # full coordinates, base 0 first
        FV = np.abs(Vf @ Fm.T)
        ok = (FV <= h[None, :]).all(1) & (np.abs(V).sum(1) > 0)
        yield FV[ok]
# cutting plane
cons = []
rng = np.random.default_rng(0)
for FV in chunks():
    if len(FV): cons.append(FV[rng.choice(len(FV), min(20, len(FV)), replace=False)])
Cn = np.vstack(cons)
for it in range(30):
    res = linprog(h.astype(float), A_ub=-Cn.astype(float), b_ub=-np.ones(len(Cn)), bounds=[(0, None)]*len(h), method='highs')
    y = res.x
    worst = []; minval = np.inf
    for FV in chunks():
        s = FV @ y
        bad = FV[s < 1 - 1e-9]
        if len(bad): worst.append(bad[np.argsort(bad @ y)[:50]])
        if len(s): minval = min(minval, s.min())
    print('iter %d: LP value %.6f, #constraints %d, min_v sum y|f(v)| = %.6f' % (it, res.fun, len(Cn), minval), flush=True)
    if not worst: break
    Cn = np.vstack([Cn] + worst)
# rationalise y and verify exactly
yr = [Fr(float(v)).limit_denominator(1000) for v in y]
den = 1
for v in yr: den = den * v.denominator // np.gcd(den, v.denominator)
yi = np.array([int(v * den) for v in yr], dtype=np.int64)
minval = min((FV @ yi).min() for FV in chunks() if len(FV))
bound = Fr(int((yi * h).sum()), den)
print('rational y (common denominator %d): min_v sum y|f(v)| = %s ; certified zeta <= %s = %.4f'
      % (den, Fr(int(minval), den), bound / Fr(int(minval), den), float(bound / Fr(int(minval), den))), flush=True)
