"""EXACT: in M8 = Syl(8,2) (T=1, m=8), dim_Q L_mixed(2)_6 < 70 <= dim_Q J(A)_6.
L_mixed(2) = J_w(R=2) + pullbacks of all families of <= 2 directions from columns u coordinates (in R = S/(x^3)).
70 <= dim_Q J(A)_6: the mod-p dimension of the (exact) line-triple local ideal is 70 (colb_syl8_out.txt), a lower bound."""
import itertools, time, numpy as np
from flint import fmpz_mat
from locrun import *
t0 = time.time()
d, cols, Ls = models.sylvester(3, 2); A = models.box(d, cols, Ls); m = 6
ring = Ring(d, m)
coords = [[1 if y == x else 0 for y in range(d)] for x in range(d)]
Dset = cols + coords
gens = width_gens(ring, width_dirs(A, d, 2, m), mod=None)
fams = [list(F) for F in itertools.combinations(Dset, 2) if np.linalg.matrix_rank(np.array(F)) == 2]
gf = family_gens(ring, A, fams, m, mod=None)
G = merge(gens, gf)
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
print(f"M8: exact dim_Q L_mixed(2)_6 = {rk} (rows {len(rows)}, N_6 = {ring.N[6]}); dim_Q J(A)_6 >= 70  => non-local at b=2: {rk < 70}  [{time.time()-t0:.0f}s]")
