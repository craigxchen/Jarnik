"""Syl(16,2), degree 7, incremental and memory-lean.  Base ideal: width powers + lines (degrees <= 7).
Since L_6 = J(A)_6 already (exact), only the DEGREE-7 generators of further families can raise L_7.
For each non-plane b-column family, project its degree-7 pullbacks to the 15 Singer components and merge into the
per-component rref basis.  Reports dim L_7 after all such families, then after adding the planes."""
import itertools, time, sys, numpy as np
from loclib import PR, image, pullback
from bigring import BigRing
from cyclib import CycRing, ideal_dims_cyc, rref_rows
import cyclib
from jfast2 import jproj_gens_fast2
from syl16 import syl16_singer
import models
t0 = time.time(); m = 7; JA7 = 9769
d, cols, Ls, colA = syl16_singer(2)
A = models.box(d, cols, Ls); ring = BigRing(d, m); cr = CycRing(ring, 15)
idx = {a: i for i, a in enumerate(colA)}
colset = {tuple(c): i for i, c in enumerate(cols)}
perm = [colset[tuple(np.roll(c, 1))] for c in cols]
gens = {}
def addpow(v):
    s = len(set((A @ np.array(v)).tolist()))
    if s <= m:
        g = ring.one()
        for t in range(s): g = ring.mul_lin(g, t, v)
        if g.any(): gens.setdefault(s, []).append(g)
addpow(cols[0]); addpow([1] * 15)
lines = sorted({tuple(sorted((a, b, a ^ b))) for a in range(1, 16) for b in range(1, 16) if a != b})
planes = [[a for a in range(1, 16) if bin(a & c).count('1') % 2 == 0] for c in range(1, 16)]
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
lreps = orbit_reps(lines); preps = orbit_reps(planes)
for ln in lreps:
    addpow(list((np.array(cols[idx[ln[0]]]) + cols[idx[ln[1]]] + cols[idx[ln[2]]]) // 2))
for V in preps:
    addpow(list(sum(np.array(cols[idx[a]]) for a in V) // 4))
for ln in lreps:
    F = [cols[idx[a]] for a in ln]
    jp, rho, fails = jproj_gens_fast2(image(A, F), 3, m)
    for i, forms in jp.items():
        for f in forms:
            v = pullback(ring, F, f, i, PR, {})
            if v.any(): gens.setdefault(i, []).append(v)
G = {k: np.array(v) for k, v in gens.items()}
orig = cyclib.rref_rows
store = []
def capture(M):
    B, rk = orig(M); store.append(B); return B, rk
cyclib.rref_rows = capture
dims = ideal_dims_cyc(cr, G, m)
cyclib.rref_rows = orig
basis = store[-15:]          # the 15 component bases in degree 7
print(f"base (width + lines): {dims}  sum of comp bases = {sum(B.shape[0] for B in basis)} [{time.time()-t0:.0f}s]", flush=True)
def merge(newvecs):
    if not len(newvecs): return
    P = cr.project(np.array(newvecs) % PR, 7)
    for s in range(15):
        M = np.vstack([basis[s], P[s]]) % PR
        basis[s], _ = rref_rows(M)
def total(): return sum(B.shape[0] for B in basis)
planeset = {tuple(sorted(idx[a] for a in V)) for V in planes}
bs = [int(x) for x in sys.argv[1].split(',')] if len(sys.argv) > 1 else [6]
for b in bs:
    seen = set(); nfam = 0
    for S in itertools.combinations(range(15), b):
        if S in seen: continue
        cur = S
        for t in range(15):
            seen.add(tuple(sorted(cur))); cur = tuple(perm[i] for i in cur)
        if any(set(S) <= set(pl) for pl in planeset): continue
        F = [cols[i] for i in S]
        P = image(A, F)
        if len(P) == np.prod([len(set(p[z] for p in P)) for z in range(b)]): continue
        jp, rho, fails = jproj_gens_fast2(P, b, m)
        cache = {}
        vecs = [pullback(ring, F, f, 7, PR, cache) for f in jp.get(7, [])]
        del cache
        merge([v for v in vecs if v.any()]); nfam += 1
        if nfam % 25 == 0: print(f"   b={b}: {nfam} families, dim L_7 = {total()} [{time.time()-t0:.0f}s]", flush=True)
    print(f"width + lines + all non-plane {b}-column families (cumulative): dim L_7 = {total()}  (J(A)_7 = {JA7}) [{time.time()-t0:.0f}s]", flush=True)
for V in preps:
    F = [cols[idx[a]] for a in V]
    jp, rho, fails = jproj_gens_fast2(image(A, F), 7, m)
    cache = {}
    vecs = [pullback(ring, F, f, 7, PR, cache) for f in jp.get(7, [])]
    del cache
    merge([v for v in vecs if v.any()])
print(f"  + planes: dim L_7 = {total()}  (J(A)_7 = {JA7}) [{time.time()-t0:.0f}s]", flush=True)
