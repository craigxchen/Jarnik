import sys, time, itertools, numpy as np, pickle
from loclib import Ring, PR, image, pullback
from cyclib import CycRing, ideal_dims_cyc
from cyc_run import dir_orbit_reps, image_sizes, powgen, fam_orbit_reps
from jfast2 import jproj_gens_fast2
import models
m = 8; t0 = time.time(); bs = [int(x) for x in sys.argv[1].split(',')]
d, cols, Ls = models.paley12(2); q = 11
A = models.box(d, cols, Ls); ring = Ring(d, m); cr = CycRing(ring, q)
JA = [0, 0, 0, 0, 0, 12, 132, 1517, 10010]
V = dir_orbit_reps(d, 1); sz = image_sizes(A, V, -22, 22)
gens = {}
for v, s in zip(V, sz):
    if s <= m:
        g = powgen(ring, v, int(s))
        if g.any(): gens.setdefault(int(s), []).append(g)
Y = A @ np.array(cols).T
cum = lambda a: list(itertools.accumulate(a))
for b in bs:
    for S in fam_orbit_reps(cols, b, q):
        P = sorted(set(map(tuple, Y[:, list(S)].tolist())))
        if len(P) == 5 ** b: continue
        F = [cols[x] for x in S]
        jp, rho, fails = jproj_gens_fast2(P, b, m)
        ng = [(i, len(f)) for i, f in jp.items() if f]
        cache = {}
        for i, forms in jp.items():
            for f in forms:
                vec = pullback(ring, F, f, i, PR, cache)
                if vec.any(): gens.setdefault(i, []).append(vec)
        print(f"  b={b} {S} |P|={len(P)} rho={rho} newgens={ng} fails={fails} [{time.time()-t0:.0f}s]", flush=True)
    G = {k: np.array(v) for k, v in gens.items()}
    Jl = ideal_dims_cyc(cr, G, m)
    print(f"== J_w + all <= {b}-col families (b in {bs[:bs.index(b)+1]}): {Jl}  J(A)={JA} cum_eq={[x==y for x,y in zip(cum(Jl),cum(JA))]} [{time.time()-t0:.0f}s]", flush=True)
