"""Numerical tangent-cone landscape of a T = 1 model (A, m), with gap (separation) constraints.  EVIDENCE ONLY.

kappa2(u) = ||(I - P_{m-1}) y_u||_2 / d^m,  y_u = (<k,u>^m)_{k in A},  P_{m-1} = orthogonal projection on Poly_{<m}|_A.
(nearstrata Prop 2.1: d^m kappa(u) is the best UNIFORM approximation error; hence
 kappa2/sqrt|A| <= kappa <= kappa2, kappa = sup_{c in V(A,m)} |Phi_c(u)|/(|c|_1 d^m).)
For each separation level s:  kappa*(s) ~ min { kappa2(u) : |u|_inf = 1 (some u_j = 1), min_{x<y} |u_x - u_y| >= s },
pairs over all n points (base point u = 0 included).  Random multistart SLSQP; reports the best value and point.
Usage: tcnum.py <profile> L m [j-list] [starts] [s-list]
"""
import sys, time, itertools, json
import numpy as np
import scipy.linalg as sla
from scipy.optimize import minimize
from flib import *
from model import profile, graded_monomials, evalmat

name = sys.argv[1]; L = int(sys.argv[2]); m = int(sys.argv[3])
jlist = [int(x) for x in sys.argv[4].split(',')] if len(sys.argv) > 4 and sys.argv[4] != 'all' else None
starts = int(sys.argv[5]) if len(sys.argv) > 5 else 6
slist = [float(x) for x in sys.argv[6].split(',')] if len(sys.argv) > 6 else [0.0, 0.02, 0.05, 0.1, 0.15, 0.2]
M, pats = profile(name)
d = M - 1
if jlist is None:
    jlist = list(range(d))
A = box(M, pats, L)
N = len(A)
t0 = time.time()
E = evalmat(A, graded_monomials(d, m - 1)).astype(float)
Qf, R, piv = sla.qr(E, mode='economic', pivoting=True)
diag = np.abs(np.diag(R))
rank = int((diag > 1e-8 * diag[0]).sum())
Q = Qf[:, :rank]
Af = np.array(A, float)
print('%s L=%d m=%d |A|=%d dim Poly_<m|A = %d, dim V(A,m) = %d (float rank; %.1fs)' % (name, L, m, N, rank, N - rank, time.time() - t0), flush=True)
dm = float(d) ** m


def resid(u):
    y = (Af @ u) ** m
    return y - Q @ (Q.T @ y)


def f(u):
    r = resid(u)
    return float(r @ r) / dm ** 2


def gf(u):
    r = resid(u)
    return 2 * m * ((r * (Af @ u) ** (m - 1)) @ Af) / dm ** 2


iu = np.triu_indices(d + 1, 1)


def seps(u):
    z = np.concatenate([[0.0], u])
    return z[iu[0]] - z[iu[1]]


def sep(u):
    return np.abs(seps(u)).min()


Jz = np.zeros((len(iu[0]), d + 1))
Jz[np.arange(len(iu[0])), iu[0]] = 1
Jz[np.arange(len(iu[0])), iu[1]] -= 1
Jz = Jz[:, 1:]

rng = np.random.default_rng(12345)
results = {}
for s in slist:
    best = (np.inf, None)
    for j in jlist:
        for st in range(starts):
            u0 = rng.uniform(-1, 1, d)
            u0[j] = 1.0
            cons = [dict(type='eq', fun=lambda u, j=j: np.array([u[j] - 1.0]), jac=lambda u, j=j: np.eye(d)[j][None, :])]
            if s > 0:
                cons.append(dict(type='ineq', fun=lambda u, s=s: seps(u) ** 2 - s ** 2,
                                 jac=lambda u: 2 * seps(u)[:, None] * Jz))
            # scale objective to O(1)
            sc = 1.0 / max(f(u0), 1e-300)
            res = minimize(lambda u: f(u) * sc, u0, jac=lambda u: gf(u) * sc, constraints=cons,
                           bounds=[(-1, 1)] * d, method='SLSQP', options=dict(maxiter=400, ftol=1e-16))
            u = res.x
            if sep(u) >= s - 1e-7 and abs(np.abs(u).max() - 1) < 1e-7:
                val = f(u) ** 0.5
                if val < best[0]:
                    best = (val, u.copy())
    results[s] = best
    u = best[1]
    if u is not None:
        z = np.sort(np.concatenate([[0.0], u]))
        print('s>=%.3f: min kappa2 = %.3e  (kappa in [%.2e, %.2e])  sep=%.4f  sorted positions=%s (%.0fs)' % (
            s, best[0], best[0] / N ** 0.5, best[0], sep(u), np.round(z, 4).tolist(), time.time() - t0), flush=True)
    else:
        print('s>=%.3f: no feasible point found' % s, flush=True)
np.save('tcnum_%s_L%d_m%d.npy' % (name, L, m), np.array([results[s][1] if results[s][1] is not None else np.full(d, np.nan) for s in slist]))
