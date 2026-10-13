"""Referee: Remark 4.2 ties in stored sign matrices. For S, decide whether some w >= 1 gives mu_xy(w) > 0 on all
D=0 pairs (LP: max t s.t. A w >= t, w >= 1, sum w = 2r). t* <= 0 <=> tie (identical cancellation)."""
import json, glob
import numpy as np
from scipy.optimize import linprog
R2 = '/tmp/claude-0/-home-user-Jarnik/ba4929a4-8665-5790-985a-3f09130a14b9/scratchpad/round7/p3data_checks/r2'
DATA = '/tmp/claude-0/-home-user-Jarnik/ba4929a4-8665-5790-985a-3f09130a14b9/scratchpad/round7/p3data_checks/data'

def tstar(S, D):
    M, r = S.shape
    P = [(x, y) for x in range(M) for y in range(x+1, M) if D[x, y] <= 0]
    if not P: return None, 0
    A = np.array([-0.5*S[x]*S[y] for x, y in P])
    c = np.zeros(r+1); c[-1] = -1
    res = linprog(c, A_ub=np.hstack([-A, np.ones((len(P), 1))]), b_ub=np.zeros(len(P)),
                  A_eq=np.hstack([np.ones((1, r)), np.zeros((1, 1))]), b_eq=[2*r], bounds=[(1, None)]*r + [(None, None)], method='highs')
    return -res.fun, len(P)

for fn in ['best_ub_M*.jsonl', 'ils2_M*.jsonl']:
    for f in sorted(glob.glob(f'{R2}/{fn}')):
        for line in open(f):
            R = json.loads(line)
            M, X, t = R['M'], R['X'], R['t']
            D = np.load(f'{DATA}/D_construct_M{M}_t{t}_u1_X{X//M}.npy')
            outs = []
            for key in ['S', 'S_fin']:
                if key in R:
                    S = np.array(R[key])
                    ts, n0 = tstar(S, D)
                    outs.append(f'{key}: r={S.shape[1]} D0pairs={n0} t*={ts:.4g}')
            print(f.split('/')[-1], M, X, t, ' | '.join(outs), flush=True)
