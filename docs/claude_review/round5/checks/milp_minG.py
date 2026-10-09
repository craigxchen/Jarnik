"""MILP (HiGHS via scipy) for m(n) = min{ G(c) : c integral, sum c = 0, ||c||_1 = n, |c_x| <= K }.
Evidence only (floating-point MILP); every reported minimiser is re-evaluated exactly."""
import sys, time
import numpy as np
from scipy.optimize import milp, LinearConstraint, Bounds
from lib import *

def minG(S, n, K=3, time_limit=60):
    M, r = S.shape
    # variables: cp[M], cm[M], z[M], t[r]
    nv = 3*M + r
    cost = np.zeros(nv); cost[3*M:] = 1.0
    A = []; lo = []; hi = []
    # t_j - S_j^T (cp - cm) >= 0 ; t_j + S_j^T(cp - cm) >= 0
    for j in range(r):
        row = np.zeros(nv); row[:M] = -S[:, j]; row[M:2*M] = S[:, j]; row[3*M+j] = 1
        A.append(row); lo.append(0); hi.append(np.inf)
        row = np.zeros(nv); row[:M] = S[:, j]; row[M:2*M] = -S[:, j]; row[3*M+j] = 1
        A.append(row); lo.append(0); hi.append(np.inf)
    row = np.zeros(nv); row[:M] = 1; row[M:2*M] = -1; A.append(row); lo.append(0); hi.append(0)
    row = np.zeros(nv); row[:2*M] = 1; A.append(row); lo.append(n); hi.append(n)
    for x in range(M):
        row = np.zeros(nv); row[x] = 1; row[2*M+x] = -K; A.append(row); lo.append(-np.inf); hi.append(0)
        row = np.zeros(nv); row[M+x] = 1; row[2*M+x] = K; A.append(row); lo.append(-np.inf); hi.append(K)
    # symmetry breaking: c -> -c ; impose z_0 ... skip
    integrality = np.zeros(nv); integrality[:3*M] = 1
    ub = np.concatenate([np.full(2*M, K), np.ones(M), np.full(r, np.inf)])
    res = milp(cost, constraints=LinearConstraint(np.array(A), lo, hi), integrality=integrality,
               bounds=Bounds(np.zeros(nv), ub), options={"time_limit": time_limit, "disp": False})
    if res.x is None:
        return None, None, res.status
    c = np.rint(res.x[:M]).astype(np.int64) - np.rint(res.x[M:2*M]).astype(np.int64)
    return Gval(S, c), c, (res.status, getattr(res, 'mip_dual_bound', None))

if __name__ == "__main__":
    q = int(sys.argv[1]); b = int(sys.argv[2]); seed = int(sys.argv[3]) if len(sys.argv) > 3 else 0
    H = normalize(paley1(q)); M = H.shape[0]
    rng = np.random.default_rng(seed)
    assign = np.concatenate([np.arange(1, M), [rng.integers(1, M)]]); rng.shuffle(assign)
    S, lab, frow = flipped_profile(H, b, assign)
    print("Paley M=%d b=%d r=%d assign=%s" % (M, b, S.shape[1], list(assign)))
    for n in range(2, 2*M + 1, 2):
        t0 = time.time()
        g, c, st = minG(S, n, K=3, time_limit=120)
        print("n=%3d  minG=%s  ratio=%.3f  F=%s  status=%s  supp=%d  maxabs=%d  (%.1fs)" % (
            n, g, g/n, F_of(H, c), st, int((c != 0).sum()), int(np.abs(c).max()), time.time()-t0), flush=True)
