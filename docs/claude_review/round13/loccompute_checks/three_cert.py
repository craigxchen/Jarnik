"""EXACT: three-cut block S1={0,1,4}, S2={0,2,3,4}, S3={1,2,3,4}, L=(1,1,1), d=5:
dim_Q L_col(3)_4 (J_w(R=2) + pullbacks of all column families incl. the full 3-column projection) = 15 < 18 <= dim_Q J(A)_4.
(18: mod-p dim of the exact mixed (b=3) local ideal, a lower bound for dim_Q J(A)_4.)"""
import itertools, numpy as np
from flint import fmpz_mat
from locrun import *
d, cols, Ls = models.cuts(5, [[0,1,4],[0,2,3,4],[1,2,3,4]], [1,1,1]); A = models.box(d, cols, Ls); m = 4
ring = Ring(d, m)
G = merge(width_gens(ring, width_dirs(A, d, 2, m), mod=None),
          *[family_gens(ring, A, [list(F) for F in itertools.combinations(cols, b)], m, mod=None) for b in (2, 3)])
rows = []
for e, gl in G.items():
    if e > m: continue
    for g in gl:
        for mu in ring.mons[m - e]:
            h = g; t = e
            for l in range(d):
                for _ in range(int(mu[l])):
                    h = ring.mul_var(h, t, l, None); t += 1
            if any(int(z) for z in h): rows.append([int(z) for z in h])
rk = fmpz_mat(rows).rank()
coords = [[1 if y == x else 0 for y in range(d)] for x in range(d)]
Dset = cols + coords
fams = [[Dset[i] for i in F] for F in itertools.combinations(range(len(Dset)), 3) if np.linalg.matrix_rank(np.array([Dset[i] for i in F])) == 3]
gm = merge(width_gens(ring, width_dirs(A, d, 2, m)), family_gens(ring, A, [list(F) for F in itertools.combinations(Dset, 2)], m), family_gens(ring, A, fams, m))
Jm = ideal_dims_modp(ring, gm, m)
print(f"exact dim_Q L_col(3)_4 = {rk};  mod-p dim of mixed(3) local ideal in degree 4 = {Jm[4]} (lower bound for dim_Q J(A)_4); non-local for columns: {rk < Jm[4]}")
