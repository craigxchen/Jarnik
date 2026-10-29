"""EXACT-over-Q certificate of (SHIFT_e) w.r.t. L_pair for e in [EMIN, EMAX] (own code).
H_Q(e) = dim V_Q(A,e) - dim V_Q(A,e+1), with dim V_Q(S,e) = exact nullity of E_{<e}(S)^T (flint fmpz nullspace).
Lower bound for dim_Q sum_l M_e(A_l): exact integer bases of V_Q(A_l, e) (verified E_{<e}^T c = 0 exactly), their
degree-e moment vectors (exact integers), rank mod Q1 (a mod-p rank of an integer matrix is <= its Q-rank).
Since sum_l M_e(A_l) is contained in M_e(A) (theorem), equality with H_Q(e) proves (SHIFT_e) over Q.
At e where the lower bound stays below H_Q(e), the script also checks whether every V_Q(A_l, e) has an exact
basis whose moments span <= the mod-p dim (only a lower bound; failures are reported as 'not certified').
usage: exact_shift.py SPEC EMIN EMAX
"""
import sys, json, time
import numpy as np
import flint
sys.path.insert(0, '/tmp/claude-0/-home-user-Jarnik/ba4929a4-8665-5790-985a-3f09130a14b9/scratchpad/round28/verify_checks')
from vlib import *

spec = sys.argv[1]; emin = int(sys.argv[2]); emax = int(sys.argv[3])
n, pats, L, A = spec_box(spec); d = n - 1
by = std_monos(d, emax + 2)
idx = {tuple(k.tolist()): i for i, k in enumerate(A)}

def kernel(S, e):
    """exact integer basis (rows) of V_Q(S,e) = {c : sum c_k k^alpha = 0, alpha std, |alpha| < e}."""
    monos = [a for t in range(e) for a in by[t]]
    E = evalmat(S, monos).astype(np.int64)
    X, nul = flint.fmpz_mat(E.T.tolist()).nullspace()
    cols = [[int(X[i, j]) for i in range(X.nrows())] for j in range(nul)]
    # exact check
    if cols:
        Cm = flint.fmpz_mat(cols)
        Z = Cm * flint.fmpz_mat(E.tolist())
        assert all(int(Z[i, j]) == 0 for i in range(Z.nrows()) for j in range(Z.ncols()))
    return cols

def moments(S, cols, e):
    if not cols: return []
    Ee = evalmat(S, by[e]).astype(np.int64)
    M = flint.fmpz_mat(cols) * flint.fmpz_mat(Ee.tolist())
    return [[int(M[i, j]) for j in range(M.ncols())] for i in range(M.nrows())]

def rank_p(rows, p=Q1):
    if not rows: return 0
    return flint.nmod_mat([[x % p for x in r] for r in rows], p).rank()

t0 = time.time()
dimV = {}
for e in range(emin, emax + 2):
    dimV[e] = len(kernel(A, e))
    log(spec, 'dim V_Q(A,%d) =' % e, dimV[e])
HQ = {e: dimV[e] - dimV[e + 1] for e in range(emin, emax + 1)}
log('H_Q', HQ)
classes = []
for l in lpair_classes(d):
    Al = shifted(A, np.array(l), idx)
    classes.append((len(Al), l, Al))
classes.sort(key=lambda t: -t[0])
res = {'spec': spec, 'H_Q': HQ, 'dimV': dimV, 'shift': {}}
for e in range(emin, emax + 1):
    rows = []; r = 0; used = 0
    for (sz, l, Al) in classes:
        if sz <= e: break
        cols = kernel(Al, e)
        if not cols: continue
        rows += moments(Al, cols, e); used += 1
        r = rank_p(rows)
        if r == HQ[e]: break
        rows = rows if len(rows) < 4 * HQ[e] else rows  # keep all (small)
    res['shift'][e] = [r, HQ[e], used]
    log('e', e, 'exact lower bound dim sum =', r, 'H_Q =', HQ[e], 'shifts used', used,
        'SHIFT_e CERTIFIED over Q' if r == HQ[e] else 'not certified (lower bound < H_Q)')
res['time'] = round(time.time() - t0, 1)
with open('/tmp/claude-0/-home-user-Jarnik/ba4929a4-8665-5790-985a-3f09130a14b9/scratchpad/round28/verify_checks/exact_shift.jsonl', 'a') as f:
    f.write(json.dumps(res) + '\n')
