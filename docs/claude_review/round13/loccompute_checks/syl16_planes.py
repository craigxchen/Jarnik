"""Syl(16,2), degrees <= 7, memory-lean rerun (the first run was OOM-killed after the plane J(P) step).
J(A) dims taken from syl16_L2_m7_out.txt (hilbert_lean, mod-p upper bounds): [0,0,0,0,0,15,540,9769].
Generators only up to the Singer shift: lines = orbit reps, planes = 1 orbit rep (the 15 planes form one orbit)."""
import itertools, time, numpy as np
from loclib import PR, image, pullback
from bigring import BigRing
from cyclib import CycRing, ideal_dims_cyc
from jfast2 import jproj_gens_fast2
from syl16 import syl16_singer
import models
t0 = time.time(); m = 7; JA = [0, 0, 0, 0, 0, 15, 540, 9769]
d, cols, Ls, colA = syl16_singer(2)
A = models.box(d, cols, Ls); ring = BigRing(d, m); cr = CycRing(ring, 15)
idx = {a: i for i, a in enumerate(colA)}
colset = {tuple(c): i for i, c in enumerate(cols)}
perm = [colset[tuple(np.roll(c, 1))] for c in cols]           # column permutation induced by the shift
gens = {}
def addpow(v):
    s = len(set((A @ np.array(v)).tolist()))
    if s <= m:
        g = ring.one()
        for t in range(s): g = ring.mul_lin(g, t, v)
        if g.any(): gens.setdefault(s, []).append(g)
addpow(cols[0]); addpow([1] * 15)
lines = sorted({tuple(sorted((a, b, a ^ b))) for a in range(1, 16) for b in range(1, 16) if a != b})
def orbit_reps(fams):
    seen = set(); reps = []
    for F in fams:
        key = tuple(sorted(idx[a] for a in F))
        if key in seen: continue
        cur = key
        for t in range(15):
            seen.add(tuple(sorted(cur))); cur = tuple(perm[i] for i in cur)
        reps.append(F)
    return reps
lreps = orbit_reps(lines)
for ln in lreps:
    z = (np.array(cols[idx[ln[0]]]) + cols[idx[ln[1]]] + cols[idx[ln[2]]]) // 2; addpow(list(z))
planes = [[a for a in range(1, 16) if bin(a & c).count('1') % 2 == 0] for c in range(1, 16)]
preps = orbit_reps(planes)
for V in preps:
    z = sum(np.array(cols[idx[a]]) for a in V) // 4; addpow(list(z))
print(f"line orbit reps {len(lreps)}, plane orbit reps {len(preps)} [{time.time()-t0:.0f}s]", flush=True)
cum = lambda a: list(itertools.accumulate(a))
def report(tag):
    G = {k: np.array(v) for k, v in gens.items()}
    Jl = ideal_dims_cyc(cr, G, m)
    print(f"  {tag}: {Jl}  J(A)={JA} cum_eq={[x==y for x,y in zip(cum(Jl),cum(JA))]} [{time.time()-t0:.0f}s]", flush=True)
report("width powers")
for ln in lreps:
    F = [cols[idx[a]] for a in ln]
    jp, rho, fails = jproj_gens_fast2(image(A, F), 3, m)
    for i, forms in jp.items():
        for f in forms:
            v = pullback(ring, F, f, i, PR, {})
            if v.any(): gens.setdefault(i, []).append(v)
report("+ lines")
for V in preps:
    F = [cols[idx[a]] for a in V]
    jp, rho, fails = jproj_gens_fast2(image(A, F), 7, m)
    cache = {}
    for i, forms in jp.items():
        for f in forms:
            v = pullback(ring, F, f, i, PR, cache)
            if v.any(): gens.setdefault(i, []).append(v)
    del cache
report("+ planes")
