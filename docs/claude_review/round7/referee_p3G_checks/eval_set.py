"""Evaluate a prime set: k<=2 closure with weight<=2 relations, then with weight-3 relations added.
usage: python3 eval_set.py primes.json M [triples]
"""
import sys, json, time
from rcore import gauss, Lprime, tight_list, gmul, conj, Ds, sigma_exact, sin2d
from rclosure import closure, status, check

def triple_list(gp, L):
    out = []; r = len(gp)
    for i in range(r):
        for j in range(i + 1, r):
            for t1 in (1, -1):
                Aij = gmul(gp[i], gp[j] if t1 > 0 else conj(gp[j]))
                for k in range(j + 1, r):
                    for t2 in (1, -1):
                        A = gmul(Aij, gp[k] if t2 > 0 else conj(gp[k]))
                        for s, D in enumerate(Ds(A)):
                            if D != 0 and L % abs(D) == 0:
                                out.append((((i, 1), (j, t1), (k, t2)), s, sigma_exact(A, s), sin2d(A, s)))
    return out

if __name__ == '__main__':
    ps = json.load(open(sys.argv[1])); M = int(sys.argv[2]); X = 2 * M
    gp = [gauss(p) for p in ps]; r = len(gp); L = Lprime(X, set(ps))
    rels = tight_list(gp, L)
    Q, log = closure(r, rels); ok, ex, hf = check(gp, L, rels, log)
    npair = sum(1 for v, *_ in rels if len(v) == 2)
    res = dict(file=sys.argv[1], r=r, tight_w2=len(rels), pair_degree=round(2 * npair / r, 2),
               k2_free=r - len(hf), k2_certs=ok)
    if len(sys.argv) > 3:
        t = time.time(); tr = triple_list(gp, L)
        rels3 = rels + tr
        Q3, log3 = closure(r, rels3); ok3, ex3, hf3 = check(gp, L, rels3, log3)
        res.update(tight_triples=len(tr), triples_per_col=round(3 * len(tr) / r, 2), k2w3_free=r - len(hf3),
                   k2w3_certs=ok3, secs=round(time.time() - t, 1))
    print(json.dumps(res), flush=True)
