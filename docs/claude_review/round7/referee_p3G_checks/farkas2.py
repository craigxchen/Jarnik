"""Arc-aware Farkas test (evidence-level numerics; exact rounding analysed separately).

Claim used: if a real z satisfies  sigma_i <v_i, z> >= 2 d_i - pi/2  for every tight relation i of a
family F, then every combination sum n_i sigma_i v_i = 0 (n_i >= 0 integers, not all 0) has
sum n_i (pi - 2 d_i) = sum n_i (pi - 2 d_i + sigma_i <v_i,z>) >= (pi/2) sum n_i >= pi
(sum n_i >= 2 since no single v_i or 2 v_i is 0). Hence no Lemma-2.4 certificate relative to {0}:
Pi_1 = {0} and the closure built from F is {0}, for ALL k and multiplicities.
LP: maximise t s.t. sigma_i <v_i,z> >= 2 d_i - pi/2 + t * [i not in an exact tie], t <= 1.
Exact ties: (v,s) and (v,s+2) both tight (then d + d' = pi/2 identically).
usage: python3 farkas2.py primes.json|random:M:b:K:seed M [w3]
"""
import sys, json, math
import numpy as np
from scipy.optimize import linprog
from scipy.sparse import lil_matrix
from rcore import gauss, Lprime, tight_list
from eval_set import triple_list
from farkas import load

def run(ps, M, w3):
    gp = [gauss(p) for p in ps]; r = len(gp); L = Lprime(2 * M, set(ps))
    rels = tight_list(gp, L) + (triple_list(gp, L) if w3 else [])
    tiedkeys = {}
    for i, (v, s, sg, s2) in enumerate(rels): tiedkeys.setdefault((v, s % 2), []).append(i)
    tie = set()
    for k, idx in tiedkeys.items():
        if len(idx) >= 2: tie.update(idx)
    n = len(rels); A = lil_matrix((n, r + 1)); b = np.zeros(n)
    for i, (v, s, sg, s2) in enumerate(rels):
        for col, co in v: A[i, col] += -sg * co
        if i not in tie: A[i, r] = 1.0
        b[i] = -(2 * math.asin(math.sqrt(float(s2))) - math.pi / 2)
    c = np.zeros(r + 1); c[r] = -1
    res = linprog(c, A_ub=A.tocsr(), b_ub=b, bounds=[(None, None)] * r + [(None, 1)], method='highs')
    out = dict(r=r, tight=n, ties=len(tie), w3=bool(w3), status=res.status)
    if res.status == 0:
        z = res.x[:r]; t = res.x[r]
        slack = A.tocsr()[:, :r] @ z - b
        out.update(t_opt=round(float(t), 6), min_viol=float(np.max(slack)) if n else 0.0, zmax=float(np.max(np.abs(z))))
    return out

if __name__ == '__main__':
    ps = load(sys.argv[1]); M = int(sys.argv[2]); w3 = len(sys.argv) > 3
    print(json.dumps(dict(set=sys.argv[1], **run(ps, M, w3))), flush=True)
