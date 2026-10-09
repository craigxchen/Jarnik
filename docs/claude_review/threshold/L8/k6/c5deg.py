"""Degree of the genuine C5-ansatz hypersurface V in (rho, sigma) space, along random lines.
Cleared determinant: row i multiplied by prod_{T ni i} A_{iT}."""
import sys, random
sys.path.insert(0, '..')
from fractions import Fraction as Fr
from polylib import *
from c5ansatz import LABELS, EDGES, TRIPLES, nbr_edges, triples_of

def cleared_matrix(rho, sigma):
    M = []
    for i in LABELS:
        e1, e2 = nbr_edges(i)
        Ts = triples_of(i)
        A = {T: (sigma[T] - rho[e1]) * (sigma[T] - rho[e2]) for T in Ts}
        row = [Fr(0)] * 5
        for k in range(3):
            T1, T2, T3 = Ts[k], Ts[(k + 1) % 3], Ts[(k + 2) % 3]
            row[TRIPLES.index(T1)] = (sigma[T2] - sigma[T3]) * A[T2] * A[T3]
        M.append(row)
    return M

def interp(xs, ys):
    n = len(xs); res = []
    for i in range(n):
        num = [Fr(1)]; den = Fr(1)
        for j in range(n):
            if j != i:
                num = mul(num, linpoly(1, -xs[j])); den *= xs[i] - xs[j]
        res = add(res, scal(ys[i] / den, num))
    return res

VARS = [('r', e) for e in EDGES] + [('s', T) for T in TRIPLES]
def point(p, v, t):
    rho = {e: p[k] + t * v[k] for k, e in enumerate(EDGES)}
    sigma = {T: p[5 + k] + t * v[5 + k] for k, T in enumerate(TRIPLES)}
    return rho, sigma

def along_line(p, v, N=30):
    xs = [Fr(k) for k in range(N)]
    ys = [det_frac(cleared_matrix(*point(p, v, x))) for x in xs]
    return interp(xs, ys)

if __name__ == '__main__':
    random.seed(int(sys.argv[1]) if len(sys.argv) > 1 else 1)
    for trial in range(4):
        p = [Fr(random.randint(-50, 50), random.randint(1, 7)) for _ in range(10)]
        v = [Fr(random.randint(-50, 50), random.randint(1, 7)) for _ in range(10)]
        D = along_line(p, v)
        rr = rational_roots(D)
        q = D; mults = {}
        for r in rr:
            while q and ev(q, r) == 0:
                q, _ = divmod_(q, linpoly(1, -r)); mults[r] = mults.get(r, 0) + 1
        import numpy as np
        print('trial', trial, 'deg D', deg(D), 'rational roots (with mult):', sum(mults.values()), 'remaining deg', deg(q))
        print('   remaining roots ~', np.round(np.sort_complex(np.roots([float(c) for c in reversed(q)])), 4))
