import sys, time, itertools, numpy as np
from loclib import Ring, PR, image, pullback
from cyclib import CycRing, ideal_dims_cyc
from cyc_run import dir_orbit_reps, image_sizes, powgen
from jfast2 import jproj_gens_fast2
import models
m = 8; t0 = time.time(); bmax = int(sys.argv[1])
d, cols, Ls = models.paley12(2); q = 11
A = models.box(d, cols, Ls); ring = Ring(d, m); cr = CycRing(ring, q)
JA = [0, 0, 0, 0, 0, 12, 132, 1517, 10010]
V = dir_orbit_reps(d, 1); sz = image_sizes(A, V, -22, 22)
gens = {}
for v, s in zip(V, sz):
    if s <= m:
        g = powgen(ring, v, int(s))
        if g.any(): gens.setdefault(int(s), []).append(g)
coords = [[1 if y == x else 0 for y in range(d)] for x in range(d)]
Dset = cols + coords
def shift(i, t): return (i // 11) * 11 + (i % 11 + t) % 11
cum = lambda a: list(itertools.accumulate(a))
Av = np.array(A)
for b in range(2, bmax + 1):
    seen = set(); stats = {}
    for S in itertools.combinations(range(22), b):
        if S in seen: continue
        for t in range(11): seen.add(tuple(sorted(shift(i, t) for i in S)))
        F = [Dset[i] for i in S]
        if np.linalg.matrix_rank(np.array(F)) < b: continue
        P = sorted(set(map(tuple, (Av @ np.array(F).T).tolist())))
        # full product of the coordinate ranges?
        full = np.prod([len(set(p[z] for p in P)) for z in range(b)])
        key = (sum(1 for i in S if i < 11), len(P) == full)
        stats[key] = stats.get(key, 0) + 1
        if len(P) == full: continue
        jp, rho, fails = jproj_gens_fast2(P, b, m)
        cache = {}
        for i, forms in jp.items():
            for f in forms:
                vec = pullback(ring, F, f, i, PR, cache)
                if vec.any(): gens.setdefault(i, []).append(vec)
    print(f"  b={b}: family orbits by (#columns, image is full product) -> {stats} [{time.time()-t0:.0f}s]", flush=True)
    G = {k: np.array(v) for k, v in gens.items()}
    Jl = ideal_dims_cyc(cr, G, m)
    print(f"== mixed (columns u coordinates) <= {b}: {Jl} J(A)={JA} cum_eq={[x==y for x,y in zip(cum(Jl),cum(JA))]} [{time.time()-t0:.0f}s]", flush=True)
