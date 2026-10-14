"""Trend survey: random twin-free profiles, d=7, L in {1,2}, P columns (random subsets).
For each: rho, ratio rho/(2L), first degree where J_w fails, b0 for columns-only and for columns+coordinates (bmax 4)."""
import sys, itertools, time, json, random
import numpy as np
from locrun import *

def b0_search(A, d, Dset, ring, m, JA, gw, bmax):
    gall = gw
    Jl = ideal_dims_modp(ring, gall, m)
    if all(cum_eq(Jl, JA)): return 1, None
    first = [i for i, e in enumerate(cum_eq(Jl, JA)) if not e][0]
    for b in range(2, bmax + 1):
        fams = [[Dset[i] for i in F] for F in itertools.combinations(range(len(Dset)), b)
                if np.linalg.matrix_rank(np.array([Dset[i] for i in F])) == b]
        gall = merge(gall, family_gens(ring, A, fams, m))
        Jl = ideal_dims_modp(ring, gall, m)
        if all(cum_eq(Jl, JA)): return b, first
    return None, first

if __name__ == '__main__':
    d = int(sys.argv[1]); ntr = int(sys.argv[2]); seed = int(sys.argv[3])
    rng = random.Random(seed)
    out = open(f'random_trend_d{d}_s{seed}_out.txt', 'w')
    for tr in range(ntr):
        L = rng.choice([1, 1, 2])
        P = rng.randint(3, 8)
        sets = []
        while len(sets) < P:
            S = tuple(sorted(rng.sample(range(d), rng.randint(2*L+1, min(d, 2*L+4)))))
            if S not in sets: sets.append(S)
        _, cols, Ls = models.cuts(d, [list(S) for S in sets], [L]*P)
        A = models.box(d, cols, Ls)
        ring0 = Ring(d, 2*d); H = hilbert_modp(A, ring0, 2*d)
        rho = max(i for i, h in enumerate(H) if h > 0); m = rho
        ring = Ring(d, m); JA = [ring.N[i]-H[i] for i in range(m+1)]
        gw = width_gens(ring, width_dirs(A, d, 2, m))
        bc, first = b0_search(A, d, [list(c) for c in cols], ring, m, JA, gw, min(4, P))
        coords = [[1 if y == x else 0 for y in range(d)] for x in range(d)]
        bm, _ = b0_search(A, d, [list(c) for c in cols] + coords, ring, m, JA, gw, 3) if bc != 1 else (1, None)
        rec = dict(L=L, P=P, sets=sets, A=len(A), rho=rho, ratio=rho/(2*L), Jw_firstfail=first, b0_col=bc, b0_mixed_le3=bm)
        line = json.dumps(rec); print(line, flush=True); out.write(line + "\n"); out.flush()
