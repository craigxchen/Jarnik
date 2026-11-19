"""Vectorised version of fullgcd.extra for k = 3."""
import numpy as np, itertools, sys
import descent_check as dc
from fullgcd import hasse_ord
K = dc.K; F27 = dc.F27; I = dc.I
def extra3(forms, signs, L, subsets=None):
    k = len(forms[0])
    vs = list(itertools.product(F27, repeat=k))
    X = np.array(vs)
    IX = K.mul(np.full(X.shape, I), X)          # i*x_b in F_729
    E = np.array([[1 + e*dc.chi27[dc.lam_val(l, x)] for x in vs] for l, e in zip(forms, signs)])
    res = {}
    subsets = subsets or [tuple(range(len(forms)))]
    fibs = []
    for alpha in range(729):
        aL = K.pow(alpha, L) if alpha else 0
        if alpha:
            val = IX[:, 0]
            for b in range(1, k):
                val = K.add[val, K.mul(IX[:, b], np.full(len(vs), K.pow(alpha, b)))]
        else:
            val = IX[:, 0]
        fib = np.nonzero(val == aL)[0]
        if len(fib) < 2: continue
        ords = np.array([hasse_ord(alpha, list(IX[idx]), L) for idx in fib])
        fibs.append((alpha, fib, ords))
    for S in subsets:
        ES = E[list(S)]; gmin = ES.min(axis=0)
        tot = 0
        for alpha, fib, ords in fibs:
            full = min(int((ES[j, fib]*ords).sum()) for j in range(len(S)))
            formal = int((gmin[fib]*ords).sum())
            tot += full - formal
        res[S] = (int(gmin.sum()), tot)
    return res
if __name__ == '__main__':
    name = sys.argv[1]; Ls = [int(a) for a in sys.argv[2:]]
    from five import fams
