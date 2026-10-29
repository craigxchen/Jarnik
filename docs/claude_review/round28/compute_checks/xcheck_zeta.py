"""Cross-check of a DFS-exhaustive 'no zonotope of order t' claim by an independent method: HiGHS MILP feasibility
restricted to the dual-reduced candidate set (sum x_v (c_v - 1) <= B - t is implied for every zonotope of order >= t).
Usage: xcheck_zeta.py SPEC t"""
import sys, math, time
sys.path.insert(0, '.')
import numpy as np
from scipy.optimize import milp, LinearConstraint, Bounds
from scipy.sparse import csr_matrix
import slib, zono
from run_box import parse
spec = sys.argv[1]; t = int(sys.argv[2])
n, pats, L = parse(spec); d = n - 1
slabs = [(S, L) for S in pats]
for cls in ('pair', 'all'):
    V, lams, bs = zono.candidates(n, slabs, cls)
    B, f, yy = zono.lp_bound(V, lams, bs)
    num, D = zono.dual_costs(V, lams, yy)
    slackN = (B - t) * D
    keep = [i for i in range(len(V)) if (int(num[i]) - D) <= slackN]
    Vk = V[keep]; ex = np.array([float(int(num[i]) - D) / D for i in keep])
    Lv = Vk @ lams.T; pos = np.maximum(Lv, 0).T; neg = np.minimum(Lv, 0).T
    nv = len(Vk)
    A1 = np.hstack([pos, lams]); A2 = np.hstack([-neg, -lams])
    rows = [A1, A2, np.concatenate([-np.ones(nv), np.zeros(d)])[None, :], np.concatenate([ex, np.zeros(d)])[None, :]]
    Acon = np.vstack(rows).astype(float)
    ub = np.concatenate([bs, bs, [-t], [float(B - t) + 1e-9]]).astype(float)
    xub = np.where(np.max(np.abs(Vk), axis=1) >= 2, 1, 2).astype(float)
    res = milp(np.zeros(nv + d), constraints=LinearConstraint(csr_matrix(Acon), -np.inf, ub),
               integrality=np.ones(nv + d), bounds=Bounds(np.concatenate([np.zeros(nv), -np.ones(d)]),
               np.concatenate([xub, np.ones(d)])), options={'time_limit': 1800})
    print(spec, cls, 'LP', B, 'reduced candidates', nv, 'MILP feasibility of order >=', t, ':', res.status, res.message, flush=True)
