"""EXACT-over-Q certificate that reach(A, L_pair) = rho (counting argument; own code).
Facts used (theorems): (i) sum_l M_e(A_l) is contained in M_e(A); (ii) sum_e H_Q(e) = |A|; (iii) a mod-p rank of an
integer matrix is <= its Q-rank.
Certificate: for every 1 <= e <= rho-1, exact integer vectors c in V_Q(U, e), U a subset of some A_l (verified
E_{<e}(U)^T c = 0 exactly), whose degree-e moment vectors have mod-p rank s_e; and dim V_Q(A, rho) = h_top,
dim V_Q(A, rho+1) = 0 (exact nullities).  If s_e = H_p(e) for e < rho and h_top = H_p(rho) and
1 + sum_{1<=e<rho} s_e + h_top = |A|, then H_Q(e) >= dim_Q sum_l M_e(A_l) >= s_e and sum H_Q = |A| force
H_Q = H_p and (SHIFT_e) over Q for all e < rho, i.e. reach_Q = rho.
usage: exact_count.py SPEC
"""
import sys, json, time, random
import numpy as np
import flint
sys.path.insert(0, '/tmp/claude-0/-home-user-Jarnik/ba4929a4-8665-5790-985a-3f09130a14b9/scratchpad/round28/verify_checks')
from vlib import *

spec = sys.argv[1]
p = Q1
n, pats, L, A = spec_box(spec); d = n - 1
t0 = time.time()
by = std_monos(d, 2 * d)
blk = least_blocks(A, by, p, 2 * d)
Hp = [b.shape[0] for b in blk]
while Hp[-1] == 0: Hp.pop()
rho = len(Hp) - 1
del blk
log(spec, '|A|', len(A), 'H_p', Hp, 'rho', rho)
idx = {tuple(k.tolist()): i for i, k in enumerate(A)}

def kernel(S, e):
    monos = [a for t in range(e) for a in by[t]]
    E = evalmat(S, monos).astype(np.int64)
    X, nul = flint.fmpz_mat(E.T.tolist()).nullspace()
    if nul == 0: return None
    Cm = flint.fmpz_mat([[X[i, j] for i in range(X.nrows())] for j in range(nul)])
    Z = Cm * flint.fmpz_mat(E.tolist())
    assert all(Z[i, j] == 0 for i in range(Z.nrows()) for j in range(Z.ncols())), 'kernel check failed'
    return Cm

def exact_nullity(S, e):
    monos = [a for t in range(e) for a in by[t]]
    E = evalmat(S, monos).astype(np.int64)
    X, nul = flint.fmpz_mat(E.T.tolist()).nullspace()
    return nul

classes = []
rng = random.Random(5)
for l in lpair_classes(d):
    Al = shifted(A, np.array(l), idx)
    perm = list(range(len(Al))); rng.shuffle(perm)
    classes.append((len(Al), l, Al[perm]))
classes.sort(key=lambda t: -t[0])
s = {}
for e in range(1, rho):
    Nlt = sum(len(by[t]) for t in range(e))
    basis = None; r = 0; used = 0
    for (sz, l, Al) in classes:
        if sz <= e: break
        m = min(sz, Nlt + Hp[e] + 200)
        U = Al[:m]
        Cm = kernel(U, e)
        if Cm is None: continue
        Ee = flint.fmpz_mat(evalmat(U, by[e]).astype(np.int64).tolist())
        Mo = Cm * Ee
        rows = [[int(Mo[i, j]) % p for j in range(Mo.ncols())] for i in range(Mo.nrows())]
        if basis is not None: rows = basis + rows
        R, r = flint.nmod_mat(rows, p).rref()
        basis = [[int(R[i, j]) for j in range(R.ncols())] for i in range(r)]
        used += 1
        if r == Hp[e]: break
    s[e] = r
    log('e', e, 'exact-vector lower bound', r, 'H_p', Hp[e], 'classes used', used)
htop = exact_nullity(A, rho); hnext = exact_nullity(A, rho + 1)
log('exact dim V_Q(A,rho) =', htop, ' dim V_Q(A,rho+1) =', hnext)
ok = all(s[e] == Hp[e] for e in range(1, rho)) and htop == Hp[rho] and hnext == 0 and \
    1 + sum(s.values()) + htop == len(A)
log('CERTIFIED reach_Q = rho =', rho if ok else None, 'ok', ok, 'time', round(time.time() - t0, 1))
with open('/tmp/claude-0/-home-user-Jarnik/ba4929a4-8665-5790-985a-3f09130a14b9/scratchpad/round28/verify_checks/exact_count.jsonl', 'a') as f:
    f.write(json.dumps({'spec': spec, 'Hp': Hp, 'rho': rho, 's': s, 'htop': htop, 'hnext': hnext, 'certified': ok,
                        'time': round(time.time() - t0, 1)}) + '\n')
