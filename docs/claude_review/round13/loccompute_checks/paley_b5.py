import sys, time, itertools, numpy as np
from loclib import Ring, PR, image, pullback
from cyclib import CycRing, ideal_dims_cyc, hilbert_cyc
from cyc_run import dir_orbit_reps, image_sizes, powgen, fam_orbit_reps
from jfast import jproj_gens_fast
import models
m = 8; t0 = time.time()
d, cols, Ls = models.paley12(2); q = 11
A = models.box(d, cols, Ls); ring = Ring(d, m); cr = CycRing(ring, q)
JA = [0, 0, 0, 0, 0, 12, 132, 1517, 10010]   # from cyc_paley_L2_out.txt (mod-p upper bounds on dim_Q J(A))
V = dir_orbit_reps(d, 1); sz = image_sizes(A, V, -22, 22)
gens = {}
for v, s in zip(V, sz):
    if s <= m:
        g = powgen(ring, v, int(s))
        if g.any(): gens.setdefault(int(s), []).append(g)
cache = {}
Y = A @ np.array(cols).T
for S in fam_orbit_reps(cols, 5, q):
    P = sorted(set(map(tuple, Y[:, list(S)].tolist())))
    if len(P) > 2000: continue
    F = [cols[x] for x in S]
    jp, rho, fails = jproj_gens_fast(P, 5, m)
    print(f"  family {S} |P|={len(P)} rho={rho} newgens {[(i, len(f)) for i, f in jp.items() if f]} fails={fails} [{time.time()-t0:.0f}s]", flush=True)
    for i, forms in jp.items():
        for f in forms:
            vec = pullback(ring, F, f, i, PR, cache)
            if vec.any(): gens.setdefault(i, []).append(vec)
G = {k: np.array(v) for k, v in gens.items()}
Jl = ideal_dims_cyc(cr, G, m, verbose=True)
cum = lambda a: list(itertools.accumulate(a))
print(f"  J_w + mod-3 congruence 5-col families: {Jl}  J(A)={JA} cum_eq={[x==y for x,y in zip(cum(Jl),cum(JA))]} [{time.time()-t0:.0f}s]")
