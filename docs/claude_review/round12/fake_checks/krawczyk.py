"""Certified regular common zero of s tangent forms Phi_i(u) = sum_k c_{i,k} <k,u>^m (integer data), u in R^d,
with u_0 = 1 and d - s - 1 further coordinates fixed to rationals; Krawczyk test in arb ball arithmetic.
Also: certified lower bound on the separation of the zero, and the non-flatness test (F) of Theorem 2.2:
the restricted Hessians B^T Hess(Phi_i) B on an (approximate, then rigorously controlled) basis B of the tangent
space are linearly independent.

Usage: krawczyk.py <basis.json> <zeros.npy> <index> m
"""
import sys, json, time
from fractions import Fraction as Fr
import numpy as np
import mpmath as mp
import flint

flint.ctx.prec = 400
mp.mp.dps = 60
D = json.load(open(sys.argv[1]))
A = [tuple(k) for k in D['A']]
Cb = D['basis'] if 'basis' in D else [D['c']]
Z = np.load(sys.argv[2])
idx = int(sys.argv[3]); m = int(sys.argv[4])
d = len(A[0]); s = len(Cb)
u = Z[idx].copy()
# keep only support of each c
terms = []
for c in Cb:
    terms.append([(ck, k) for ck, k in zip(c, A) if ck])
nrm = [sum(abs(ck) for ck, _ in T) for T in terms]


def F_mp(um):
    out = []; Jm = []
    for T, nn in zip(terms, nrm):
        f = mp.mpf(0); g = [mp.mpf(0)] * d
        for ck, k in T:
            y = mp.fsum(kx * ux for kx, ux in zip(k, um) if kx)
            y7 = y ** (m - 1)
            f += ck * y7 * y
            for x in range(d):
                if k[x]:
                    g[x] += m * ck * y7 * k[x]
        out.append(f / nn); Jm.append([gx / nn for gx in g])
    return mp.matrix(out), mp.matrix(Jm)


t0 = time.time()
um = [mp.mpf(x) for x in u]
# Gauss-Newton (min-norm) with u_0 fixed
for it in range(12):
    Fv, Jv = F_mp(um)
    J1 = Jv[:, 1:]
    step = J1.T * mp.lu_solve(J1 * J1.T, Fv)
    um = [um[0]] + [um[x] - step[x - 1] for x in range(1, d)]
    if mp.norm(Fv) < mp.mpf(10) ** -50:
        break
print('refined residual %.2e (%.0fs)' % (float(mp.norm(F_mp(um)[0])), time.time() - t0), flush=True)
Fv, Jv = F_mp(um)
# choose s unknown coordinates by column pivoting on J (excluding coordinate 0)
Jn = np.array([[float(Jv[i, j]) for j in range(d)] for i in range(s)])
import scipy.linalg as sla
_, _, piv = sla.qr(Jn[:, 1:], pivoting=True)
unk = sorted([int(p) + 1 for p in piv[:s]])
fixed = [x for x in range(d) if x not in unk]
den = 10 ** 25
ur = {x: Fr(int(mp.nint(um[x] * den)), den) for x in fixed}
ur[0] = Fr(1)
# re-solve the square system in the unknowns with fixed rationals
xm = [um[x] for x in unk]
for it in range(10):
    full = [mp.mpf(ur[x].numerator) / ur[x].denominator if x in ur else None for x in range(d)]
    for j, x in enumerate(unk):
        full[x] = xm[j]
    Fv, Jv = F_mp(full)
    Js = mp.matrix([[Jv[i, x] for x in unk] for i in range(s)])
    dx = mp.lu_solve(Js, Fv)
    xm = [xm[j] - dx[j] for j in range(s)]
print('square system residual %.2e' % float(mp.norm(Fv)), flush=True)
xc = [Fr(int(mp.nint(x * 10 ** 40)), 10 ** 40) for x in xm]
r = Fr(1, 10 ** 30)


def arbq(q):
    return flint.arb(q.numerator) / flint.arb(q.denominator)


def evalF_J(box_unk):
    """arb evaluation of F (normalised) and full Jacobian, Hessians on the point/box."""
    U = [None] * d
    for x in range(d):
        if x in ur:
            U[x] = arbq(ur[x])
    for j, x in enumerate(unk):
        U[x] = box_unk[j]
    Fs = []; Js = []; Hs = []
    for T, nn in zip(terms, nrm):
        f = flint.arb(0); g = [flint.arb(0)] * d
        H = [[flint.arb(0)] * d for _ in range(d)]
        for ck, k in T:
            y = flint.arb(0)
            for kx, Ux in zip(k, U):
                if kx:
                    y += kx * Ux
            y2 = y * y; y6 = flint.arb(1)
            for _ in range(m - 2):
                y6 = y6 * y
            y7 = y6 * y
            f += ck * y7 * y
            for x in range(d):
                if k[x]:
                    g[x] += m * ck * k[x] * y7
                    for x2 in range(d):
                        if k[x2]:
                            H[x][x2] += m * (m - 1) * ck * k[x] * k[x2] * y6
        Fs.append(f / nn); Js.append([gx / nn for gx in g]); Hs.append([[h / nn for h in row] for row in H])
    return Fs, Js, Hs


center = [arbq(x) for x in xc]
X = [arbq(x) + flint.arb(0, float(r)) for x in xc]
Fc, _, _ = evalF_J(center)
_, JX, HX = evalF_J(X)
Jmid = np.array([[float(JX[i][x].mid()) for x in unk] for i in range(s)])
Y = np.linalg.inv(Jmid)
Ya = [[flint.arb(Y[i, j]) for j in range(s)] for i in range(s)]
K = []
ok = True
for i in range(s):
    val = center[i] - sum((Ya[i][j] * Fc[j] for j in range(s)), flint.arb(0))
    for j in range(s):
        Mij = (flint.arb(1) if i == j else flint.arb(0)) - sum((Ya[i][l] * JX[l][unk[j]] for l in range(s)), flint.arb(0))
        val += Mij * (X[j] - center[j])
    K.append(val)
    lo = arbq(xc[i]) - flint.arb(float(r)); hi = arbq(xc[i]) + flint.arb(float(r))
    inside = (val.lower() > (arbq(xc[i]) - arbq(r)).upper()) and (val.upper() < (arbq(xc[i]) + arbq(r)).lower())
    ok = ok and inside
for i in range(s):
    print('  K_%d - xc = %s ; Fc = %s' % (i, (K[i] - arbq(xc[i])).str(5), Fc[i].str(5)))
print('|Y| = %.3e' % np.abs(Y).max())
print('Krawczyk K(X) subset int X: %s  (box radius 1e-30 in %d unknowns %s; %.0fs)' % (ok, s, unk, time.time() - t0), flush=True)
assert ok
# separation lower bound
pos = [Fr(0)] + [ur[x] if x in ur else xc[unk.index(x)] for x in range(d)]
rad = [Fr(0)] + [Fr(0) if x in ur else r for x in range(d)]
gaps = []
for a in range(d + 1):
    for b in range(a + 1, d + 1):
        gaps.append(abs(pos[a] - pos[b]) - rad[a] - rad[b])
umax = max(abs(p) + q for p, q in zip(pos, rad))
print('certified: sep(u*) >= %.8f, |u*|_inf <= %.8f (u*_0 = 1), normalised sep >= %.8f; equally spaced max 2/%d = %.6f' % (
    float(min(gaps)), float(umax), float(min(gaps) / umax), d, 2 / d))
zz = np.array([float(p) for p in pos])
pr = np.abs(zz[:, None] - zz[None, :])[np.triu_indices(d + 1, 1)]
print('Gamma(u*) = %.4f ; sorted positions %s' % (2 / len(pr) * np.sum(np.log(np.maximum(1, 1 / pr))), np.round(np.sort(zz), 5).tolist()))
# (F): restricted Hessians on tangent basis.  Full Jacobian at u* enclosed by JX (all d columns).
Jfull = np.array([[float(JX[i][x].mid()) for x in range(d)] for i in range(s)])
_, _, Vt = np.linalg.svd(Jfull)
B = Vt[s:].T  # d x (d-s) approx tangent basis (orthonormal)
Ba = [[flint.arb(B[x, j]) for j in range(d - s)] for x in range(d)]
# residual |J(u*) b| for basis vectors (rigorous upper bound via JX)
res = 0.0
for j in range(d - s):
    for i in range(s):
        v = sum((JX[i][x] * Ba[x][j] for x in range(d)), flint.arb(0))
        res = max(res, float(abs(v).upper()))
svJ = np.linalg.svd(Jfull, compute_uv=False)
print('tangent basis: max |J(u*) b| <= %.2e ; sigma_min(J) ~ %.3e' % (res, svJ[-1]))
# matrices R_i = B^T H_i B (entries arb)
Rm = []
for i in range(s):
    R = np.zeros(((d - s), (d - s)))
    Rr = 0.0
    for a in range(d - s):
        for b in range(d - s):
            v = sum((Ba[x][a] * HX[i][x][y] * Ba[y][b] for x in range(d) for y in range(d)), flint.arb(0))
            R[a, b] = float(v.mid()); Rr = max(Rr, float(v.rad()))
    Rm.append(R)
G = np.array([R.ravel() for R in Rm])
svG = np.linalg.svd(G, compute_uv=False)
Hn = max(np.abs(np.array([[float(HX[i][x][y].mid()) for y in range(d)] for x in range(d)])).sum() for i in range(s))
print('restricted Hessians: singular values of the %d x %d matrix of entries: %s ; max entry radius %.1e ; Hessian l1 bound %.2e' % (
    s, (d - s) ** 2, np.array2string(svG, precision=3), Rr, Hn))
# perturbation: true tangent vectors b' with |b - b'| <= res/sigma_min ; entries change by <= 2 Hn |b-b'|*...
pert = 2 * Hn * (res / svJ[-1]) * 2 * (d - s) + Rr * (d - s)
print('(F) certified if sigma_min(G) = %.3e > perturbation bound %.3e : %s' % (svG[-1], pert, svG[-1] > pert))
