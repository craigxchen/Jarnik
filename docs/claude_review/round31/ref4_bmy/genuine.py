"""Run the lattice-level BMY engine and (IV.2) on genuine function-field clusters (record examples)."""
import sys, itertools
import sympy as sp
from fractions import Fraction as Fr
from bmy_global import analyse, t

def cluster(roots, K):
    D = len(roots)
    Q = sp.expand(sp.prod([t-r for r in roots]))
    pts = []
    for F in itertools.product([0, 1], repeat=D):
        z = sp.expand(sp.prod([(t-sp.conjugate(roots[k])) if F[k] else (t-roots[k]) for k in range(D)]))
        dz = sp.expand(z-Q)
        if dz == 0 or sp.degree(dz, t) <= K:
            pts.append(F)
    # verify pairwise
    zs = [sp.expand(sp.prod([(t-sp.conjugate(roots[k])) if F[k] else (t-roots[k]) for k in range(D)])) for F in pts]
    pw = max(sp.degree(sp.expand(a-b), t) for a, b in itertools.combinations(zs, 2))
    assert pw <= K
    cols, nl = [], []
    for F in pts:
        P = sp.expand(sp.prod([t-roots[k] for k in range(D) if not F[k]]))
        re, im = sp.expand((P+sp.conjugate(P))/2), sp.expand((P-sp.conjugate(P))/(2*sp.I))
        # conjugate of polynomial with real t: sympy conjugate treats t as complex; t declared real below
        cols.append((re, im)); nl.append(sp.expand(P*sp.conjugate(P)))
    N = sp.expand(Q*sp.conjugate(Q))
    return cols, N, nl

if __name__ == '__main__':
    x0 = sp.Rational(1, 3)
    for y, K in [([1, 2, 3, 4, 5, -6, -7], 4), ([1, 2, 3, 4, -5, -6, -7, 8], 5)]:
        roots = [x0+sp.I*v for v in y]
        cols, N, nl = cluster(roots, K)
        res = analyse(cols, N, nl)
        print('D=%d K=%d' % (len(y), K), {k: (str(v) if isinstance(v, Fr) else v) for k, v in res.items()})
