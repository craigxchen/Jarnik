"""Abstract Sylvester PLANE image: P = Lambda_3 cap [-L,L]^7, Lambda_3 = Z-span of the simplex-code vectors u_x
(the image of the 7 columns of a 3-dim subspace when every coset is large).  Coordinates c (y = U c), so P_c in Z^7.
Compare, in Q[c_1..c_7], degrees <= rho+1:  J(P)  vs  J_w(P) (directions v in [-R,R]^7 in c-coords)
                                              vs  J_w + pullbacks of the 7 LINE families (b=3).   Z_7 (Singer of F_8) isotypic, mod p."""
import sys, itertools, time, numpy as np
from sympy import Matrix
from flint import nmod_mat
from loclib import PR, image, pullback
from fullring import FullRing
from cyclib import CycRing, ideal_dims_cyc
from jfast2 import jproj_gens_fast2
t0 = time.time(); L = int(sys.argv[1]); R = int(sys.argv[2])
# F_8 = F_2[t]/(t^3+t+1); points alpha^l; columns a = alpha^l too (so the shift permutes both)
pts = []; x = 1
for l in range(7):
    pts.append(x); x <<= 1
    if x & 8: x ^= 0b1011
U = Matrix([[bin(a & y).count('1') % 2 for y in pts] for a in pts])   # rows: columns a, cols: points x
D = int(U.det()); Uin = np.array(U.inv() * D, dtype=np.int64)
Y = np.array(list(itertools.product(range(-L, L + 1), repeat=7)), dtype=np.int64)
C = Y @ Uin.T; ok = (C % D == 0).all(axis=1)
Pc = C[ok] // D                      # points in c-coordinates
Umat = np.array(U, dtype=np.int64)
print(f"L={L}: |P| = {len(Pc)}, det U = {D} [{time.time()-t0:.0f}s]", flush=True)
# Hilbert function of P (full ring), mod p, via pivots of transposed evaluation
Dmax = 2 * 7 * L
H = []
ring0 = FullRing(7, 14)
basis = np.zeros((0, len(Pc)), dtype=np.int64)
Pco = Pc.astype(object)
for i in range(0, 15):
    M = ring0.mons[i]
    before = basis.shape[0]
    for c0 in range(0, len(M), 1500):
        Mb = M[c0:c0 + 1500]
        V = np.ones((len(Mb), len(Pc)), dtype=object)
        for z in range(7):
            V = V * (Pco[:, z][None, :] ** Mb[:, z][:, None].astype(object))
        V = np.array([[int(v) % PR for v in r] for r in V], dtype=np.int64)
        allr = np.vstack([basis, V])
        Mt = nmod_mat(allr.shape[0], allr.shape[1], [int(z) for z in allr.ravel()], PR)
        Rr, rk = Mt.rref()
        basis = np.array([int(z) for z in Rr.entries()], dtype=np.int64).reshape(allr.shape)[:rk]
    H.append(basis.shape[0] - before)
    if basis.shape[0] == len(Pc): break
rho = len(H) - 1
print(f"  H_P = {H}  rho = {rho} [{time.time()-t0:.0f}s]", flush=True)
m = rho + 1
ring = FullRing(7, m); cr = CycRing(ring, 7)
JP = [ring.N[i] - (H[i] if i < len(H) else 0) for i in range(m + 1)]
# shift check: coordinate shift of c corresponds to x -> alpha x; P_c invariant?
Ps = set(map(tuple, Pc.tolist()))
assert all(tuple(np.roll(p, 1)) in Ps for p in Ps), "P_c not shift invariant"
gens = {}
seen = set(); nd = 0
for v in itertools.product(range(-R, R + 1), repeat=7):
    if not any(v) or v in seen: continue
    orb = set()
    for t in range(7):
        orb.add(tuple(np.roll(v, t))); orb.add(tuple(-z for z in np.roll(v, t)))
    seen |= orb
    import math
    if math.gcd(*[abs(z) for z in v if z]) != 1: continue
    s = len(set((Pc @ np.array(v)).tolist()))
    if s <= m:
        g = ring.one()
        for t in range(s): g = ring.mul_lin(g, t, v)
        if g.any(): gens.setdefault(s, []).append(g); nd += 1
cumf = lambda a: list(itertools.accumulate(a))
G = {k: np.array(vv) for k, vv in gens.items()}
Jw = ideal_dims_cyc(cr, G, m)
print(f"  J(P)        = {JP}\n  J_w(R={R}) [{nd} orbit reps] = {Jw}  cum_eq={[a==b for a,b in zip(cumf(Jw), cumf(JP))]} [{time.time()-t0:.0f}s]", flush=True)
# lines: triples of columns {a, b, a^b}; in c-coords the column forms are the rows of U
lines = sorted({tuple(sorted((a, b, a ^ b))) for a in pts for b in pts if a != b})
pos = {a: i for i, a in enumerate(pts)}
for ln in lines:
    F = [list(Umat[pos[a]]) for a in ln]
    Pimg = image(Pc, F)
    jp, rh, fails = jproj_gens_fast2(Pimg, 3, m)
    for i, forms in jp.items():
        for f in forms:
            vec = pullback(ring, F, f, i, PR, {})
            if vec.any(): gens.setdefault(i, []).append(vec)
G = {k: np.array(vv) for k, vv in gens.items()}
Jl = ideal_dims_cyc(cr, G, m)
print(f"  J_w + 7 line pullbacks = {Jl}  cum_eq={[a==b for a,b in zip(cumf(Jl), cumf(JP))]} [{time.time()-t0:.0f}s]", flush=True)
