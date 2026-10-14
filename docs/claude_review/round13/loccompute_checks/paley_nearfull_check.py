"""For every non-grid 5-column Paley family that is NOT a congruence family (|P| > 2000): mod-P62 cumulative dims of J(P)
in degrees <= 7 equal those of the grid ideal (y_q^5).  Since dim_{P62} J(P)_{<=i} >= dim_Q J(P)_{<=i} >= dim_Q (y^5)_{<=i},
equality proves J(P)_{<=7} = (y_q^5)_{<=7} over Q (so these families add nothing in degrees <= 7)."""
import itertools, numpy as np, models
from flint import nmod_mat
from jfast import P62
from loclib import bmonos
from cyc_run import fam_orbit_reps
d, cols, Ls = models.paley12(2); A = models.box(d, cols, Ls); Y = A @ np.array(cols).T
allok = True
for S in fam_orbit_reps(cols, 5, 11):
    P = np.array(sorted(set(map(tuple, Y[:, list(S)].tolist()))))
    if len(P) < 2000: continue
    mons = [e for t in range(8) for e in bmonos(5, t)]
    E = np.ones((len(P), len(mons)), dtype=object)
    for j, e in enumerate(mons):
        for z in range(5):
            if e[z]: E[:, j] = E[:, j] * P[:, z].astype(object) ** e[z]
    rk = nmod_mat(E.shape[0], E.shape[1], [int(x) % P62 for x in E.ravel()], P62).rank()
    dimJ = len(mons) - rk
    grid = sum(1 for e in mons if max(e) >= 5)
    ok = dimJ == grid; allok &= ok
    print(S, len(P), "dim_P62 J(P)_{<=7} =", dimJ, " grid ideal dim =", grid, ok, flush=True)
print("ALL OK:", allok)
