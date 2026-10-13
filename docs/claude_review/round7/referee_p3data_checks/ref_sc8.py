"""Referee: own semi-continuous MILP W_sc at M=8 (all 127 cuts): min sum lam_s, mu >= D, lam_s in {0} u [w0, U]."""
import json, sys, math, time
import numpy as np
from scipy.optimize import milp, LinearConstraint, Bounds
DATA = '/tmp/claude-0/-home-user-Jarnik/ba4929a4-8665-5790-985a-3f09130a14b9/scratchpad/round7/p3data_checks/data'
M = 8
cuts = []
for mask in range(1, 2**(M-1)):
    s = np.array([-1 if (mask >> x) & 1 else 1 for x in range(M)]); cuts.append(s)
S = np.array(cuts)  # n x M
iu = np.triu_indices(M, 1)
A = (0.5 - 0.5*S[:, iu[0]]*S[:, iu[1]] - 0.5).T  # pairs x cuts
n = S.shape[0]
w0s = {16: math.log(17), 32: math.log(37)}
tl = float(sys.argv[1]) if len(sys.argv) > 1 else 300
out = []
for X in (16, 32):
    for t in range(6):
        D = np.load(f'{DATA}/D_construct_M8_t{t}_u1_X{X//8}.npy'); d = D[iu]
        w0 = w0s[X]; U = 4*d.max() + 4*w0
        # vars lam (n), z (n)
        c = np.concatenate([np.ones(n), np.zeros(n)])
        rows = [np.hstack([A, np.zeros((len(d), n))])]; lo = [d]; hi = [np.full(len(d), np.inf)]
        rows.append(np.hstack([np.eye(n), -U*np.eye(n)])); lo.append(np.full(n, -np.inf)); hi.append(np.zeros(n))
        rows.append(np.hstack([np.eye(n), -w0*np.eye(n)])); lo.append(np.zeros(n)); hi.append(np.full(n, np.inf))
        con = LinearConstraint(np.vstack(rows), np.concatenate(lo), np.concatenate(hi))
        integ = np.concatenate([np.zeros(n), np.ones(n)])
        t0 = time.time()
        res = milp(c, constraints=con, integrality=integ, bounds=Bounds(np.zeros(2*n), np.concatenate([np.full(n, U), np.ones(n)])),
                   options=dict(time_limit=tl, mip_rel_gap=1e-7))
        lam = res.x[:n] if res.x is not None else None
        ncut = int((lam > 1e-7).sum()) if lam is not None else None
        rec = dict(X=X, t=t, W=res.fun, dual_bound=getattr(res, 'mip_dual_bound', None), status=res.status, msg=res.message[:40], ncuts=ncut,
                   lam2=2*np.linalg.eigvalsh(D)[-1], secs=time.time()-t0)
        out.append(rec); print(json.dumps(rec), flush=True)
json.dump(out, open('ref_sc8_out.json', 'w'), indent=0)
