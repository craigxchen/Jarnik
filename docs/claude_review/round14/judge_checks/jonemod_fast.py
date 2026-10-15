"""Judge check J2 (fast variant: int32, v and -v merged) (round 14): one-modulus test at q0 for ps, K, in three classes, own conventions (jglib).
System: genuine labels (all levels) at every q in F \\ {q0}; free moduli collision-free (no demand); at q0 an
arbitrary parity-respecting level-1 label ell (level-1 relaxation, as round 10 / both referees).
Classes: PU (per-unit budget m^2 <= N*2^[s odd]); ANG (PU + exact per-v box test, see jang.py);
SLg (PU + slab system at phi = arg pi_j, which is exactly m_(v,s) <= |D_s(Z_v)| for every v, s).
Counts: all survivors / outside the genuine class (c = u L mod n) / outside every eps-class (c = u eps L mod n).
Usage: python3 jonemod.py q0 K p1 p2 ... (F read from referee forced_K*.txt, given as env FFILE).
"""
import sys, os, math, itertools
import numpy as np
from fractions import Fraction as Fr
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from jglib import two_sq, genuine, nu, Zv, Ds, VK, vq
from jang import box_meets_circle

q0 = int(sys.argv[1]); K = int(sys.argv[2]); ps = [int(x) for x in sys.argv[3:]]
F = None
for line in open(os.environ['FFILE']):
    if line.startswith('F '): F = eval(line[2:])
assert q0 in F
r = len(ps); pis = [two_sq(p) for p in ps]; V = VK(r, K)
m = q0-1 if q0 % 4 == 1 else q0+1; n = m//4
flags = {'PU': [], 'ANG': [], 'SLg': []}
for v in V:
    z = Zv(pis, v); N = z[0]**2 + z[1]**2
    base = []
    for s in range(4):
        d = Ds(z, s); g = 1
        for q in F:
            if q != q0: g *= q**vq(d, q)
        base.append(g)
    for key in flags: flags[key].append([False]*4)
    for s in range(4):  # flag = adding one factor q0 at unit s breaks the class
        ms = list(base); ms[s] *= q0
        pu = all(ms[t]**2 <= N*(2 if t % 2 else 1) for t in range(4))
        ang = pu and box_meets_circle(Fr(2*ms[2]**2, N)-1, 1-Fr(2*ms[0]**2, N), Fr(ms[3]**2, N)-1, 1-Fr(ms[1]**2, N))
        slg = pu and all(ms[t] <= abs(Ds(z, t)) for t in range(4))
        flags['PU'][-1][s] = not pu; flags['ANG'][-1][s] = not ang; flags['SLg'][-1][s] = not slg
    # base system itself must be fine in every class (sanity)
    assert all(base[t]**2 <= N*(2 if t % 2 else 1) for t in range(4))
par = nu(ps, q0)
axes = [np.arange(par[j], m, 2, dtype=np.int64) for j in range(r)]
grid = np.stack(np.meshgrid(*axes, indexing='ij'), axis=-1).reshape(-1, r).astype(np.int32)
idx = {v: i for i, v in enumerate(V)}
bad = {key: np.zeros(len(grid), dtype=bool) for key in flags}
for i, v in enumerate(V):
    if next(x for x in v if x) < 0: continue  # -v handled with v: <-v,ell> = -<v,ell>, unit s -> -s
    j = idx[tuple(-x for x in v)]
    d = (grid @ np.array(v, dtype=np.int32)) % m
    hit = (d % n == 0); s = (d // n) % 4
    for key in flags:
        fl = np.array([flags[key][i][t] or flags[key][j][(-t) % 4] for t in range(4)]); bad[key] |= hit & fl[s]
_, L = genuine(ps, q0)
units = [u for u in range(1, n+1) if math.gcd(u, n) == 1] if n > 1 else [1]
gen = set(tuple((u*x) % n for x in L) for u in units)
epsc = set(tuple((u*e*x) % n for e, x in zip(eps, L)) for u in units for eps in itertools.product([1, -1], repeat=r))
out = {}
for key in flags:
    surv = grid[~bad[key]]
    cs = [tuple(int(x) % n for x in row) for row in surv]
    ex = [tuple(int(x) for x in row) for row, c in zip(surv, cs) if c not in epsc][:3]
    out[key] = (len(surv), sum(c not in gen for c in cs), sum(c not in epsc for c in cs), ex)
print(f"ps={ps} K={K} q0={q0} |F|={len(F)} labels={len(grid)} genuine L={L}")
for key, (a, b, c, ex) in out.items():
    print(f"  {key}: survivors {a} / outside genuine class {b} / outside every eps-class {c}; examples {ex}")
