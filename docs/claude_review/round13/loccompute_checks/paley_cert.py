"""EXACT certificate (over Q) that J(A)_7 is not contained in the local ideal L_7, Paley-12, L=2, T=1, m=8.
Works in the Z_11-invariant part: P_0 = orbit-sum coordinates (rational).
 (1) rows P_0(mu*g), g exact integer generators of L (degree <= 7), mu monomials: exact fmpz rank r0 = dim_Q P_0(L_7).
 (2) exact invariant polynomials f (deg <= 7) vanishing on A (verified at ALL points of A over Z); top forms T_f.
 (3) rank_p([rows; P_0(T_f)]) > r0  ==> some T_f (in J(A)_7) is not in L_7 over Q.   [mod-p rank <= Q rank]
usage: python3 paley_cert.py [cong]   (cong: include the exact pullbacks of the six mod-3 congruence 5-column families)"""
import sys, time, itertools, math, numpy as np
from flint import nmod_mat, fmpz_mat
from loclib import Ring, PR, pullback
from cyclib import CycRing
from cyc_run import dir_orbit_reps, image_sizes, fam_orbit_reps
from jfast import P62, ratrec
from jfast2 import jproj_gens_fast2
import models
t0 = time.time(); D = 7
d, cols, Ls = models.paley12(2); q = 11
A = models.box(d, cols, Ls); ring = Ring(d, D); cr = CycRing(ring, q)
def P0(h, i):   # exact orbit-sum coordinates of an object/int vector over R_i
    O = cr.orb[i]
    return [int(sum(int(h[j]) for j in O[o])) for o in range(O.shape[0])]
def powex(v, e):
    vec = np.ones(1, dtype=object)
    for t in range(e): vec = ring.mul_lin(vec, t, v, None)
    return vec
gens = []   # (degree, exact object vector)
V = dir_orbit_reps(d, 1); sz = image_sizes(A, V, -22, 22)
for v, s in zip(V, sz):
    if s <= D:
        g = powex(v, int(s))
        if any(int(z) for z in g): gens.append((int(s), g))
print("width gens (orbit reps, exact):", [(e, sum(1 for x in gens if x[0]==e)) for e in sorted({x[0] for x in gens})], flush=True)
tag = 'J_w(R=1) = L_cols(4)'
if 'cong' in sys.argv:
    tag = 'J_w + all 5-col families (only the congruence ones are not full boxes)'
    Y = A @ np.array(cols).T
    for S in fam_orbit_reps(cols, 5, q):
        P = sorted(set(map(tuple, Y[:, list(S)].tolist())))
        if len(P) > 2000: continue          # near-full: new gens only y_q^5 (checked in paley_bgen_out.txt)
        F = [cols[x] for x in S]
        jp, rho, fails = jproj_gens_fast2(P, 5, D)
        for i, forms in jp.items():
            for f in forms:
                g = pullback(ring, F, f, i, None, {})
                if any(int(z) for z in g): gens.append((i, g))
    print("  with congruence families:", len(gens), "gens", flush=True)
# (1) rows
rows = []
for e, g in gens:
    for mu in ring.mons[D - e]:
        h = g
        t = e
        for l in range(d):
            for _ in range(int(mu[l])):
                h = ring.mul_var(h, t, l, None); t += 1
        rows.append(P0(h, D))
rows = [r for r in rows if any(r)]
r0 = fmpz_mat(rows).rank()
print(f"(1) {tag}: #rows {len(rows)}, EXACT dim_Q P_0(L_7) = {r0}  (invariant coords: {cr.orb[D].shape[0]}) [{time.time()-t0:.0f}s]", flush=True)
# (2) invariant vanishing polynomials of degree <= 7
basis = [(0, None)] + [(i, o) for i in range(1, D + 1) for o in range(cr.orb[i].shape[0])]
def monvals(M, Pts):
    vals = np.ones((len(M), len(Pts)), dtype=np.int64)
    for x in range(d):
        e = M[:, x][:, None]; kx = Pts[:, x][None, :]
        vals *= np.where(e == 0, 1, np.where(e == 1, kx, kx * kx))
    return vals
cols_all = [np.ones(len(A), dtype=np.int64)]
for i in range(1, D + 1):
    Vi = monvals(ring.mons[i], A)
    O = cr.orb[i]
    for o in range(O.shape[0]):
        cols_all.append(Vi[O[o]].sum(axis=0))
E = np.stack(cols_all, axis=1)          # (|A|, nb) exact int64
# point-orbit reps
seen = set(); reps = []
for a, k in enumerate(A.tolist()):
    k = tuple(k)
    if k in seen: continue
    for t in range(q): seen.add(tuple(np.roll(k, t)))
    reps.append(a)
Er = E[reps]
from crtker import kernel_exact, verify
from echelon import Echelon
ntop = cr.orb[D].shape[0]; nlow = E.shape[1] - ntop
Rp = [[x % PR for x in r] for r in rows]
ech = Echelon(ntop)
for r in Rp: ech.add(r)
cur = ech.rank
def select(v0, p0):
    # cheap pre-test mod p0 (top part outside the current span, tested mod PR after reduction mod p0 is NOT valid,
    # so only require a nonzero top part here; the real test is done after exact lifting)
    return any(v0[nlow:])
found = []
for iv in kernel_exact(Er, nprimes=4, select=select):
    tv = [z % PR for z in iv[nlow:]]
    if any(ech.reduce(tv)):
        ok = verify(E, iv)   # exact, ALL points of A
        if ok:
            ech.add(tv); found.append((iv, tv))
print(f"(2) invariant polys deg<=7: {E.shape[1]} unknowns; exact lifted kernel vectors with new top forms: {len(found)} [{time.time()-t0:.0f}s]", flush=True)
rk = nmod_mat(len(Rp) + len(found), ntop, [x for r in Rp + [f[1] for f in found] for x in r], PR).rank()
print(f"(3) exact vanishing invariant polys found: {len(found)}; rank_p(rows + their top forms) = {rk} > r0 = {r0}: {rk > r0}", flush=True)
if found:
    iv = found[0][0]
    print("    first certificate polynomial: max |coef| =", max(abs(z) for z in iv), " nonzero top coeffs:", sum(1 for z in iv[nlow:] if z))
    open('paley_cert_poly%s.txt' % ('_cong' if 'cong' in sys.argv else ''), 'w').write(repr(iv))
print(f"VERDICT: {'J(A)_7 not contained in L_7 over Q (rigorous)' if rk > r0 else 'no certificate'} [{time.time()-t0:.0f}s]")
