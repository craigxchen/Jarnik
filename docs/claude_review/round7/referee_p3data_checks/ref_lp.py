"""Referee: own W_cut LP (column generation, brute-force max-cut pricing in numpy) with certified dual
bound <nu,D>/surplus(nu), surplus computed by exhaustive enumeration; 2 lambda_max; compare with table42_rows.json."""
import json, sys, time
import numpy as np
from scipy.optimize import linprog
DATA = '/tmp/claude-0/-home-user-Jarnik/ba4929a4-8665-5790-985a-3f09130a14b9/scratchpad/round7/p3data_checks/data'
ROWS = '/tmp/claude-0/-home-user-Jarnik/ba4929a4-8665-5790-985a-3f09130a14b9/scratchpad/round7/p3data_checks/r2/table42_rows.json'

def cuts_chunk(M, lo, hi):
    # cuts indexed by integers 0..2^(M-1)-1 with s_{M-1} = +1 fixed; bits give s_0..s_{M-2}
    idx = np.arange(lo, hi, dtype=np.int64)
    bits = ((idx[:, None] >> np.arange(M-1)) & 1).astype(np.int8)
    S = np.ones((hi-lo, M), dtype=np.int8); S[:, :M-1] = 1 - 2*bits
    return S

def best_cuts(M, Nmat, k=30):
    """max over nonconstant cuts of sum_{x<y} nu_xy([sep]-1/2) = -(1/4) s^T N s (N symmetric zero diag, N_xy=nu_xy).
    returns surplus and top-k cuts."""
    tot = 2**(M-1); ch = 1 << 16
    best = []; bestval = -np.inf
    for lo in range(1, tot, ch):  # skip index 0 = constant cut
        hi = min(tot, lo + ch)
        S = cuts_chunk(M, lo, hi).astype(np.float64)
        val = -0.25*np.einsum('ij,ij->i', S @ Nmat, S)
        if len(val) > k:
            sel = np.argpartition(-val, k)[:k]
        else:
            sel = np.arange(len(val))
        for i in sel: best.append((val[i], lo + i))
        best.sort(key=lambda t: -t[0]); best = best[:k]
    S = cuts_chunk(M, 0, tot)[[b[1] for b in best]] if False else np.array([cuts_chunk(M, b[1], b[1]+1)[0] for b in best])
    return best[0][0], S

def wcut(D, verbose=False):
    M = D.shape[0]
    iu = np.triu_indices(M, 1); d = D[iu]
    # initial pool: random cuts
    rng = np.random.default_rng(1)
    pool = {}
    for _ in range(4*M):
        s = rng.choice([-1, 1], size=M); s[-1] = 1
        if abs(s.sum()) < M: pool[tuple(s)] = s
    it = 0
    while True:
        it += 1
        S = np.array(list(pool.values()), dtype=np.float64)  # n x M
        A = 0.5 - 0.5*(S[:, iu[0]]*S[:, iu[1]])  # [sep]
        A = (A - 0.5).T  # pairs x cuts : mu coefficient
        # min 1.lam s.t. A lam >= d  -> -A lam <= -d ; add slack-free; may be infeasible early -> add big penalty slack
        n = S.shape[0]; P = len(d)
        c = np.concatenate([np.ones(n), 1e3*np.ones(P)])
        Aub = -np.hstack([A, np.eye(P)])
        res = linprog(c, A_ub=Aub, b_ub=-d, bounds=(0, None), method='highs')
        nu = -res.ineqlin.marginals  # duals >= 0
        Nmat = np.zeros((M, M)); Nmat[iu] = nu; Nmat = Nmat + Nmat.T
        surplus, top = best_cuts(M, Nmat)
        val = res.fun
        lb = nu @ d / surplus if surplus > 0 else np.inf
        if verbose: print(it, n, val, lb, surplus, flush=True)
        if surplus <= 1 + 1e-9 or it > 400:
            lam = res.x[:n]; slack = res.x[n:]
            return val, lb, slack.max()
        for s in top:
            pool[tuple(s)] = s

if __name__ == '__main__':
    rows = json.load(open(ROWS))
    Ms = [int(a) for a in sys.argv[1:]] or [8, 12, 20]
    out = []
    for R in rows:
        M, X, t = R['M'], R['X'], R['t']
        if M not in Ms: continue
        D = np.load(f'{DATA}/D_construct_M{M}_t{t}_u1_X{X//M}.npy')
        t0 = time.time()
        lam2 = 2*np.linalg.eigvalsh(D)[-1]
        W, lb, sl = wcut(D)
        rec = dict(M=M, X=X, t=t, lam2=lam2, Wcut=W, Wlb=lb, maxslack=sl, auth_lam2=R['lam2'], auth_Wcut=R['Wcut'],
                   Dmax=D.max(), auth_Dmax=R['Dmax'], secs=time.time()-t0)
        out.append(rec)
        print(json.dumps(rec), flush=True)
    json.dump(out, open('ref_lp_out_%s.json' % '_'.join(map(str, Ms)), 'w'), indent=1)
