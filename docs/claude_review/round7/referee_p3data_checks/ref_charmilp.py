"""Referee: own MILP for min s_c = sum_j w_j |v_j(c)| - W/2 over zero-sum integer c != 0, |c_x| <= B,
v = A^T c, A = (1+S)/2, for the realised profiles; also exact min pair margin mu (over all pairs) and over D=0 pairs."""
import json, glob, math, sys
import numpy as np
from scipy.optimize import milp, LinearConstraint, Bounds
R2 = '/tmp/claude-0/-home-user-Jarnik/ba4929a4-8665-5790-985a-3f09130a14b9/scratchpad/round7/p3data_checks/r2'
DATA = '/tmp/claude-0/-home-user-Jarnik/ba4929a4-8665-5790-985a-3f09130a14b9/scratchpad/round7/p3data_checks/data'

def decode(mask, M):
    return np.array([-1 if (mask >> x) & 1 else 1 for x in range(M)])

def charmin(cols, M, B=12, tl=120):
    S = np.array([decode(m, M) for m, _ in cols]).T
    A = (1 + S)//2; r = A.shape[1]
    w = np.array([sum(math.log(p) for p in pl) for _, pl in cols]); W = w.sum()
    # vars: c (M ints), t (r cont)
    nv = M + r
    cost = np.concatenate([np.zeros(M), w])
    rows = []; lo = []; hi = []
    # t_j - v_j >= 0 ; t_j + v_j >= 0
    for j in range(r):
        a = np.zeros(nv); a[:M] = -A[:, j]; a[M+j] = 1; rows.append(a); lo.append(0); hi.append(np.inf)
        a = np.zeros(nv); a[:M] = A[:, j]; a[M+j] = 1; rows.append(a); lo.append(0); hi.append(np.inf)
    a = np.zeros(nv); a[:M] = 1; rows.append(a); lo.append(0); hi.append(0)
    Acon = np.array(rows)
    integrality = np.concatenate([np.ones(M), np.zeros(r)])
    best = (np.inf, None)
    for x0 in range(M):
        lb = np.concatenate([-B*np.ones(M), np.zeros(r)]); ub = np.concatenate([B*np.ones(M), np.full(r, np.inf)])
        lb[x0] = 1
        # symmetry: also require c_y <= 0 ... no; just c_x0>=1 (some coord positive)
        res = milp(cost, constraints=LinearConstraint(Acon, lo, hi), integrality=integrality,
                   bounds=Bounds(lb, ub), options=dict(time_limit=tl, mip_rel_gap=1e-9))
        if res.x is not None and res.fun - W/2 < best[0]:
            best = (res.fun - W/2, np.round(res.x[:M]).astype(int), res.status)
    return best, S, w

if __name__ == '__main__':
    which = sys.argv[1] if len(sys.argv) > 1 else 'all'
    out = []
    for fn, key in [('ils2_M*.jsonl', 'prime_ils_cols'), ('best_ub_M*.jsonl', 'prime_cols')]:
        for f in sorted(glob.glob(f'{R2}/{fn}')):
            for line in open(f):
                r = json.loads(line)
                if key not in r or r['M'] == 8: continue
                if which != 'all' and str(r['M']) != which: continue
                M, X, t = r['M'], r['X'], r['t']
                D = np.load(f'{DATA}/D_construct_M{M}_t{t}_u1_X{X//M}.npy')
                (val, c, st), S, w = charmin(r[key], M)
                W = w.sum()
                mu = np.array([[0.5*(w[S[x] != S[y]].sum() - w[S[x] == S[y]].sum()) for y in range(M)] for x in range(M)])
                iu = np.triu_indices(M, 1)
                minmu = mu[iu].min(); minmu0 = mu[iu][D[iu] == 0].min() if (D[iu] == 0).any() else None
                rec = dict(src=f.split('/')[-1], M=M, X=X, t=t, charmin=val, supp=int((c != 0).sum()), c=c.tolist(),
                           status=int(st), min_pair_mu=minmu, min_mu_D0=minmu0)
                out.append(rec); print(json.dumps({k: v for k, v in rec.items() if k != 'c'}), flush=True)
    json.dump(out, open(f'ref_charmilp_{which}.json', 'w'), indent=0)
