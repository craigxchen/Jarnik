"""Syl(16,2), degrees <= 6: width powers (columns, sigma, line half-sums, plane quarter-sums) + ALL 2-dim families from
columns u coordinates.  J(A) dims from syl16_L2_m7_out.txt (mod-p upper bounds): [0,0,0,0,0,15,540]."""
import itertools, time, numpy as np
from loclib import PR, image, pullback
from bigring import BigRing
from cyclib import CycRing, ideal_dims_cyc
from jfast2 import jproj_gens_fast2
from syl16 import syl16_singer
import models
t0 = time.time(); m = 6; JA = [0, 0, 0, 0, 0, 15, 540]
d, cols, Ls, colA = syl16_singer(2)
A = models.box(d, cols, Ls); ring = BigRing(d, m); cr = CycRing(ring, 15)
gens = {}
def addpow(v):
    s = len(set((A @ np.array(v)).tolist()))
    if s <= m:
        g = ring.one()
        for t in range(s): g = ring.mul_lin(g, t, v)
        if g.any(): gens.setdefault(s, []).append(g)
addpow(cols[0]); addpow([1] * 15)
idx = {a: i for i, a in enumerate(colA)}
for a in range(1, 16):
    for b in range(a + 1, 16):
        c = a ^ b
        if c > b: addpow(list((np.array(cols[idx[a]]) + cols[idx[b]] + cols[idx[c]]) // 2))
coords = [[1 if y == x else 0 for y in range(d)] for x in range(d)]
Dset = cols + coords
Av = np.array(A); ngrid = 0; nnon = 0
for S in itertools.combinations(range(30), 2):
    F = [Dset[i] for i in S]
    if np.linalg.matrix_rank(np.array(F)) < 2: continue
    P = sorted(set(map(tuple, (Av @ np.array(F).T).tolist())))
    full = np.prod([len(set(p[z] for p in P)) for z in range(2)])
    if len(P) == full: ngrid += 1; continue
    nnon += 1
    jp, rho, fails = jproj_gens_fast2(P, 2, m)
    for i, forms in jp.items():
        for f in forms:
            vec = pullback(ring, F, f, i, PR, {})
            if vec.any(): gens.setdefault(i, []).append(vec)
print(f"2-dim families: {ngrid} grids, {nnon} non-grid [{time.time()-t0:.0f}s]", flush=True)
G = {k: np.array(v) for k, v in gens.items()}
Jl = ideal_dims_cyc(cr, G, m)
print(f"width + all mixed pairs: {Jl}  J(A) = {JA}  [{time.time()-t0:.0f}s]")
