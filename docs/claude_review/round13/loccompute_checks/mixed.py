"""Minimal b for families from a direction set Dset (columns + coordinates [+ extra]) on small models."""
import sys, itertools, time, json
import numpy as np
from locrun import *

def mixed_b0(d, cols, Ls, m=None, R=2, extra=(), bmax=None, verbose=True):
    A = models.box(d, cols, Ls)
    ring0 = Ring(d, 2*d); H = hilbert_modp(A, ring0, 2*d)
    rho = max(i for i, h in enumerate(H) if h > 0)
    m = m or rho
    ring = Ring(d, m)
    JA = [ring.N[i]-H[i] for i in range(m+1)]
    coords = [[1 if y == x else 0 for y in range(d)] for x in range(d)]
    Dset = [list(c) for c in cols] + coords + [list(e) for e in extra]
    gall = width_gens(ring, width_dirs(A, d, R, m))
    out = {'|A|': len(A), 'rho': rho, 'm': m}
    for b in range(1, (bmax or len(Dset)) + 1):
        if b > 1:
            fams = []
            for F in itertools.combinations(range(len(Dset)), b):
                M = np.array([Dset[i] for i in F])
                if np.linalg.matrix_rank(M) < b: continue
                fams.append([Dset[i] for i in F])
            gall = merge(gall, family_gens(ring, A, fams, m))
        Jl = ideal_dims_modp(ring, gall, m)
        eq = cum_eq(Jl, JA)
        if verbose: print(f"   b={b}: L={Jl} J(A)={JA} eq={eq}", flush=True)
        if all(eq):
            out['b0_mixed'] = b; return out
    out['b0_mixed'] = None
    return out

if __name__ == '__main__':
    for sets, Lv in [([[0,1,4],[0,2,3,4],[1,2,3,4]], (1,1,1)), ([[0,1],[0,2,3],[1,2,3]], (1,1,2)),
                     ([[0,1,3,4,5],[0,2,3,4,5],[1,2,3,4,5]], (1,1,1))]:
        d = 1 + max(max(S) for S in sets)
        print(sets, Lv, flush=True)
        print(mixed_b0(*models.cuts(d, sets, Lv)))
