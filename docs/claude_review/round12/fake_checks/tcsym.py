"""Symmetry-reduced numerical tangent-cone landscape (EVIDENCE ONLY), same quantity as tcnum.py:
kappa2(u) = ||(I - P_{<m}) y_u||_2 / d^m, computed blockwise in the isotypic decomposition of a cyclic coordinate
permutation group <sigma> of order g preserving A (orthonormal Fourier basis on orbits).
Usage: tcsym.py <profile> L m <perm: 'shift' | 'mult:a'> starts s-list
"""
import sys, time
import numpy as np
import scipy.linalg as sla
from scipy.optimize import minimize
from flib import *
from model import profile

name = sys.argv[1]; L = int(sys.argv[2]); m = int(sys.argv[3]); pspec = sys.argv[4]
starts = int(sys.argv[5]); slist = [float(x) for x in sys.argv[6].split(',')]
M, pats = profile(name)
d = M - 1; q = d
if pspec == 'shift':
    perm = [(x + 1) % q for x in range(q)]
else:
    a = int(pspec.split(':')[1]); perm = [(a * x) % q for x in range(q)]
g = 1; cur = perm[:]
while cur != list(range(q)):
    cur = [perm[x] for x in cur]; g += 1


def act(k):
    out = [0] * d
    for x, v in enumerate(k):
        out[perm[x]] = v
    return tuple(out)


A = box(M, pats, L)
Aset = set(A); assert all(act(k) in Aset for k in A)
t0 = time.time()
# orbits of points: list of (rep, [sigma^j rep])
seen = set(); porb = []
for k in A:
    if k in seen:
        continue
    o = [k]
    for _ in range(g - 1):
        o.append(act(o[-1]))
    seen |= set(o); porb.append(o)
mons = []
for j in range(m):
    mons += reduced_monomials(d, j)
seen = set(); morb = []
for a_ in mons:
    if a_ in seen:
        continue
    o = [a_]
    for _ in range(g - 1):
        o.append(act(o[-1]))
    seen |= set(o); morb.append(a_)
Mon = np.array(morb, dtype=np.int64)
w = np.exp(2j * np.pi / g)
free = [i for i, o in enumerate(porb) if len(set(o)) == g]
fixed = [i for i, o in enumerate(porb) if len(set(o)) == 1]
assert len(free) + len(fixed) == len(porb), 'only free and fixed orbits supported (g prime)'
# index arrays: point index of sigma^j rep
pidx = {k: i for i, k in enumerate(A)}
IF = np.array([[pidx[porb[i][j]] for j in range(g)] for i in free], dtype=np.int64)
IX = np.array([pidx[porb[i][0]] for i in fixed], dtype=np.int64)
Aa = np.array(A, dtype=np.int64)


def monvals(P):
    V = np.ones((len(P), len(Mon)))
    for x in range(d):
        ax = Mon[:, x][None, :]; px = P[:, x][:, None].astype(float)
        V = V * np.where(ax == 0, 1.0, np.where(ax == 1, px, px * px))
    return V


Qs = []
for s in range(g):
    if s == 0:
        rowsF = sum(monvals(Aa[IF[:, j]]) for j in range(g)) / np.sqrt(g)
        rows = np.vstack([rowsF, monvals(Aa[IX])]) if len(IX) else rowsF
    else:
        rows = sum(monvals(Aa[IF[:, j]]) * w ** (-s * j) for j in range(g)) / np.sqrt(g)
    Q, R, P = sla.qr(rows, mode='economic', pivoting=True)
    dg = np.abs(np.diag(R)); rk = int((dg > 1e-8 * dg[0]).sum())
    Qs.append(Q[:, :rk])
    print('block s=%d: %d x %d, rank %d -> dim V_s = %d' % (s, rows.shape[0], rows.shape[1], rk, rows.shape[0] - rk), flush=True)
print('total dim V(A,%d) = %d (%.0fs)' % (m, sum(Q.shape[0] - Q.shape[1] for Q in Qs), time.time() - t0), flush=True)
Af = Aa.astype(float)
dm = float(d) ** m
Wm = np.array([[w ** (-s * j) for j in range(g)] for s in range(g)]) / np.sqrt(g)


def resid_full(u):
    y = (Af @ u) ** m
    YF = y[IF]  # free x g
    Yh = YF @ Wm.T  # free x g  (column s)
    r = np.zeros(len(A))
    tot = 0.0
    RF = np.zeros_like(Yh)
    for s in range(g):
        v = Yh[:, s]
        if s == 0 and len(IX):
            v = np.concatenate([v, y[IX]])
        rr = v - Qs[s] @ (Qs[s].conj().T @ v)
        tot += float(np.vdot(rr, rr).real)
        if s == 0 and len(IX):
            r[IX] = rr[len(IF):].real
            rr = rr[:len(IF)]
        RF[:, s] = rr
    # back to point values: r(sigma^j rep) = sum_s conj(Wm[s,j]) * RF[:, s]
    r[IF] = (RF @ Wm.conj()).real
    return tot, r


def f(u):
    return resid_full(u)[0] / dm ** 2


def gf(u):
    _, r = resid_full(u)
    return 2 * m * ((r * (Af @ u) ** (m - 1)) @ Af) / dm ** 2


# sanity: gradient check
u0 = np.random.default_rng(0).uniform(-1, 1, d)
e = np.zeros(d); e[0] = 1e-6
print('grad check: %.4e vs %.4e' % ((f(u0 + e) - f(u0 - e)) / 2e-6, gf(u0)[0]), flush=True)
iu = np.triu_indices(d + 1, 1)
Jz = np.zeros((len(iu[0]), d + 1)); Jz[np.arange(len(iu[0])), iu[0]] = 1; Jz[np.arange(len(iu[0])), iu[1]] -= 1; Jz = Jz[:, 1:]


def seps(u):
    z = np.concatenate([[0.0], u]); return z[iu[0]] - z[iu[1]]


rng = np.random.default_rng(12345)
jl = [0] if pspec == 'shift' else list(range(d))
for s in slist:
    best = (np.inf, None)
    for j in jl:
        for st in range(starts):
            u0 = rng.uniform(-1, 1, d); u0[j] = 1.0
            cons = [dict(type='eq', fun=lambda u, j=j: np.array([u[j] - 1.0]), jac=lambda u, j=j: np.eye(d)[j][None, :])]
            if s > 0:
                cons.append(dict(type='ineq', fun=lambda u, s=s: seps(u) ** 2 - s ** 2, jac=lambda u: 2 * seps(u)[:, None] * Jz))
            sc = 1.0 / max(f(u0), 1e-300)
            res = minimize(lambda u: f(u) * sc, u0, jac=lambda u: gf(u) * sc, constraints=cons, bounds=[(-1, 1)] * d,
                           method='SLSQP', options=dict(maxiter=400, ftol=1e-16))
            u = res.x
            if np.abs(seps(u)).min() >= s - 1e-7 and abs(np.abs(u).max() - 1) < 1e-7:
                val = f(u) ** 0.5
                if val < best[0]:
                    best = (val, u.copy())
    u = best[1]
    if u is None:
        print('s>=%.3f: none' % s, flush=True); continue
    z = np.sort(np.concatenate([[0.0], u]))
    print('s>=%.3f: min kappa2 = %.3e (kappa in [%.2e, %.2e]) sep=%.4f sorted=%s (%.0fs)' % (
        s, best[0], best[0] / len(A) ** 0.5, best[0], np.abs(seps(u)).min(), np.round(z, 4).tolist(), time.time() - t0), flush=True)
