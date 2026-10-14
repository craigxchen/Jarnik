"""For the six mod-3 congruence 5-column families: the ideal generated (mod PR) by the EXACT verified generators returned by
jfast2 equals, degree by degree up to 7, the mod-P62 dimension of J(P).  Sandwich: dim_PR(gen. ideal) <= dim_Q(gen. ideal)
<= dim_Q J(P) <= dim_P62 J(P) (cumulative).  Equality => the exact generators span J_Q(P)_{<=7}, so L_col(5) rows used in
paley_cert3 are complete."""
import itertools, numpy as np, models
from flint import nmod_mat
from jfast import P62
from jfast2 import jproj_gens_fast2
from loclib import bmonos, PR
from fullring import FullRing
from cyc_run import fam_orbit_reps
d, cols, Ls = models.paley12(2); A = models.box(d, cols, Ls); Y = A @ np.array(cols).T
fr = FullRing(5, 7); allok = True
for S in fam_orbit_reps(cols, 5, 11):
    P = sorted(set(map(tuple, Y[:, list(S)].tolist())))
    if len(P) > 2000: continue
    jp, rho, fails = jproj_gens_fast2(P, 5, 7)
    # generated ideal dims mod PR
    dims = []; basis = []
    for i in range(8):
        rows = [fr.mul_var(r, i - 1, l) for r in basis for l in range(5)] if i else []
        for f in jp.get(i, []):
            v = np.zeros(fr.N[i], dtype=np.int64)
            for e, c in f.items(): v[fr.idx[i][int(np.array(e) @ fr.pw)]] = c % PR
            rows.append(v)
        if not rows: dims.append(0); basis = []; continue
        M = np.array(rows) % PR
        Mt = nmod_mat(M.shape[0], M.shape[1], [int(z) for z in M.ravel()], PR); R_, rk = Mt.rref()
        basis = list(np.array([int(z) for z in R_.entries()], dtype=np.int64).reshape(M.shape)[:rk]); dims.append(rk)
    # P62 dims of J(P)
    Pa = np.array(P, dtype=object); jd = []; prev = 0
    rows = []
    for i in range(8):
        for e in bmonos(5, i):
            col = np.ones(len(P), dtype=object)
            for z in range(5):
                if e[z]: col = col * Pa[:, z] ** e[z]
            rows.append([int(x) % P62 for x in col])
        rk = nmod_mat(len(rows), len(P), [x for r in rows for x in r], P62).rank()
        jd.append(len(bmonos(5, i)) - (rk - prev)); prev = rk
    ok = list(itertools.accumulate(dims)) == list(itertools.accumulate(jd)); allok &= ok
    print(S, len(P), "generated (mod PR):", dims, " J(P) (mod P62):", jd, ok, flush=True)
print("ALL OK:", allok)
