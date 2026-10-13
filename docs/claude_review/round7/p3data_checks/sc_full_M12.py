"""Full semi-continuous MILP at M = 12 (all 2047 cuts, indicator binaries; HiGHS with a time limit):
a lower bound (MILP dual bound) and an upper bound for W_sc, to calibrate the rwl1 heuristic of q_lean.py.
The big-M constant is U = Wub, a known feasible value (q_lean.py's rwl1 bound): every solution with
W <= Wub has all cut weights <= Wub, so the MILP dual bound is a valid lower bound for W_sc.
Usage: python3 sc_full_M12.py t Xf time_limit Wub
"""
import sys, math, json
import numpy as np
from qsolve import lp_full, sc_milp
from qcore import primes_1mod4_above

t, Xf, tl = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
Wub = float(sys.argv[4])
M = 12
D = np.load('data/D_construct_M%d_t%d_u1_X%d.npy' % (M, t, Xf))
wmin = math.log(primes_1mod4_above(Xf * M, 3)[0])
Wc, lam, masks, nu = lp_full(D)
U = Wub
r = sc_milp(D, masks, wmin, U=U, time_limit=tl)
lamv = r['lam']
rec = dict(M=M, t=t, X=Xf * M, Wcut=Wc, U=U, status=int(r['status']), sc_ub=r['W'], sc_lb=r['lb'],
           ncuts=int((lamv > 0).sum()) if lamv is not None else None,
           max_weight=float(lamv.max()) if lamv is not None else None, time=r['time'])
print(json.dumps(rec))
with open('sc_full_M12.jsonl', 'a') as f:
    f.write(json.dumps(rec) + '\n')
