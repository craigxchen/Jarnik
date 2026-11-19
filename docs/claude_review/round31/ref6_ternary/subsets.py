"""For q=27 (and q=3), k=2,3: best formal common degree g over all m-subsets of the (3^k-1)/2 forms and all signs."""
import numpy as np, itertools, sys
from f9lib import Fq
def forms_list(k):
    out = []
    for c in itertools.product(range(3), repeat=k):
        if any(c) and c[[i for i in range(k) if c[i]][0]] == 1:
            out.append(c)
    return out
def chi_arrays(r, k):
    F = Fq(r); q = F.q
    vecs = np.array(list(itertools.product(range(q), repeat=k)))
    forms = forms_list(k)
    addt = np.array(F.addt); neg = np.array(F.neg); chi = np.array(F.chi)
    res = []
    for lam in forms:
        s = np.zeros(len(vecs), np.int64)
        for c, col in zip(lam, vecs.T):
            if c == 1: s = addt[s, col]
            elif c == 2: s = addt[s, neg[col]]
        res.append(chi[s])
    return forms, np.array(res), q**k, q**(k-1)
if __name__ == '__main__':
    r, k, mmax = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
    forms, C, n, h = chi_arrays(r, k)
    print('q=3^%d k=%d n=%d h=%d #forms=%d' % (r, k, n, h, len(forms)))
    for m in range(2, mmax+1):
        best = (-1,)
        seen = 0
        for S in itertools.combinations(range(len(forms)), m):
            if S[0] != 0: break   # WLOG (GL_3(F_3) transitive on forms) first form fixed
            for eps in itertools.product((1, -1), repeat=m-1):
                e = np.array((1,) + eps)[:, None]
                g = int(np.min(1 + e*C[list(S)], axis=0).sum())
                if g > best[0]: best = (g, [forms[s] for s in S], (1,)+eps)
        print(m, 'best g =', best[0], 'g-h =', best[0]-h, best[1:], flush=True)
