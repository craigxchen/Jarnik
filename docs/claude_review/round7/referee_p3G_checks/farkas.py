"""Farkas certificate that a prime set has NO Lemma-2.4 certificate relative to Pi_0 = {0}.

If an integer vector y satisfies sigma_i <v_i, y> > 0 for every tight relation (v_i, s_i) of a family F,
then no combination sum n_i sigma_i v_i (n_i >= 1, relations from F) is 0, so Lemma 2.4 / Def. 2.5
gives Pi_1 = {0}, hence Pi* = {0}, for the closure built from F -- for ANY k, n and arcs.
Also tests the converse direction for random sets (LP infeasible <=> some nonneg combination is 0).
usage: python3 farkas.py primes.json|random:M:b:K:seed M [w3]
"""
import sys, json, random
import numpy as np
from scipy.optimize import linprog
from scipy.sparse import lil_matrix
from rcore import gauss, Lprime, tight_list, rand_primes_like_author
from eval_set import triple_list

def load(spec):
    if spec.startswith('random:'):
        _, M, b, K, seed = spec.split(':'); M, b, K, seed = int(M), int(b), float(K), int(seed)
        Kv = int(K) if K == int(K) else K; P = int(M ** Kv)
        return rand_primes_like_author(b * (M - 1), P, 2 * P, random.Random(seed))
    return json.load(open(spec))

def farkas(ps, M, w3):
    gp = [gauss(p) for p in ps]; r = len(gp); L = Lprime(2 * M, set(ps))
    rels = tight_list(gp, L) + (triple_list(gp, L) if w3 else [])
    A = lil_matrix((len(rels), r))
    for i, (v, s, sg, s2) in enumerate(rels):
        for col, co in v: A[i, col] += -sg * co      # -sigma <v,y> <= -1
    res = linprog(np.zeros(r), A_ub=A.tocsr(), b_ub=-np.ones(len(rels)), bounds=[(None, None)] * r,
                  method='highs')
    out = dict(r=r, tight=len(rels), w3=bool(w3), lp_status=res.status, lp_msg=res.message[:40])
    if res.status == 0:
        y = res.x
        for scale in (1, 10, 100, 1000, 10 ** 4, 10 ** 6):
            yi = [int(round(t * scale)) for t in y]
            if all(sum(co * yi[col] for col, co in v) * sg > 0 for v, s, sg, s2 in rels):
                out.update(exact_certificate=True, scale=scale, max_abs_y=max(abs(t) for t in yi))
                break
        else:
            out.update(exact_certificate=False)
    return out

if __name__ == '__main__':
    ps = load(sys.argv[1]); M = int(sys.argv[2]); w3 = len(sys.argv) > 3
    print(json.dumps(dict(set=sys.argv[1], **farkas(ps, M, w3))), flush=True)
