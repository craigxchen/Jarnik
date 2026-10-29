"""2-D: S = S' u (S'+l0) with S' a Sidon-type set (all differences distinct). Compare shift reach, zeta, and
(mod-P) graded Hilbert function. Then exact (sympy, over Q) re-check of (SHIFT_e) for the reported range."""
import sys, itertools, random, time
sys.path.insert(0, '/tmp/claude-0/-home-user-Jarnik/ba4929a4-8665-5790-985a-3f09130a14b9/scratchpad/round27/interp_checks')
from shiftlib import *
from t4_search2d import graded_H, shift_reach, zeta
import sympy as sp
random.seed(int(sys.argv[1]) if len(sys.argv) > 1 else 3)
G = 5
cand = [(a, b) for a in range(G + 1) for b in range(G + 1)]
random.shuffle(cand)
Sp = []
diffs = set()
for p in cand:
    nd = set()
    ok = True
    for q in Sp:
        for dd in [(p[0]-q[0], p[1]-q[1]), (q[0]-p[0], q[1]-p[1])]:
            if dd in diffs or dd in nd: ok = False; break
            nd.add(dd)
        if not ok: break
    if ok:
        Sp.append(p); diffs |= nd
    if len(Sp) >= 7: break
l0 = (7, 3)
S = sorted(set(Sp) | set((a + l0[0], b + l0[1]) for a, b in Sp))
print("S' =", sorted(Sp), 'l0 =', l0, '|S| =', len(S))
print('graded H_S (mod P):', graded_H(S, 2), " graded H_S':", graded_H(Sp, 2))
ms, rho = shift_reach(S, 2)
print('shift criterion certifies ZE^all(S,m) for m <=', ms, '(rho =', rho, ')', flush=True)
z = zeta(S)
print('zeta(S) (max order of binomial-product form found in S; numeric search) =', z, flush=True)
# exact check over Q of SHIFT_e for e < ms, and of failure at e = ms (if ms <= rho)
def M_exact(Sx, e):
    low = [(i, t - i) for t in range(e) for i in range(t + 1)]
    top = [(i, e - i) for i in range(e + 1)]
    El = sp.Matrix([[k[0]**a * k[1]**b for (a, b) in low] for k in Sx]) if low else None
    Ee = sp.Matrix([[k[0]**a * k[1]**b for (a, b) in top] for k in Sx])
    if El is None:
        ker = [sp.eye(len(Sx))[:, i] for i in range(len(Sx))]
    else:
        ker = El.T.nullspace()
    if not ker: return sp.zeros(0, e + 1)
    Kmat = sp.Matrix.hstack(*ker).T
    return Kmat * Ee
Sset = set(S)
Ls = sorted(set((a[0]-b[0], a[1]-b[1]) for a in S for b in S) - {(0, 0)})
for e in range(0, min(ms, rho) + 1):
    MA = M_exact(S, e); rA = MA.rank() if MA.shape[0] else 0
    rows = []
    for l in Ls:
        Al = [k for k in S if (k[0]+l[0], k[1]+l[1]) in Sset]
        if len(Al) <= e: continue
        Ml = M_exact(Al, e)
        if Ml.shape[0]: rows.append(Ml)
    rs = sp.Matrix.vstack(*rows).rank() if rows else 0
    print('  exact: e', e, 'dim M_e(S)', rA, 'dim sum_l M_e(S_l)', rs, 'SHIFT' if rA == rs else 'FAIL', flush=True)
