"""sdc11 m=8: min over {u_0=1, |u|_inf<=1, sep>=s} of F(u)=sum_i (Phi_i(u)/(|c_i|_1 d^8))^2 (5 forms). Evidence."""
import json, time, sys, numpy as np
from scipy.optimize import minimize
D = json.load(open('sdc11_m8_basis.json'))
A = np.array(D['A'], float); C = np.array(D['basis'], float); C /= np.abs(C).sum(1, keepdims=True)
d = 11; m = 8; dm = float(d) ** m
def Ph(u): return C @ (A @ u) ** m / dm
def F(u): p = Ph(u); return float(p @ p)
def G(u):
    y = A @ u; J = m * (C * y ** (m - 1)) @ A / dm
    return 2 * Ph(u) @ J
iu = np.triu_indices(d + 1, 1)
Jz = np.zeros((len(iu[0]), d + 1)); Jz[np.arange(len(iu[0])), iu[0]] = 1; Jz[np.arange(len(iu[0])), iu[1]] -= 1; Jz = Jz[:, 1:]
def gaps(u): z = np.concatenate([[0.0], u]); return z[iu[0]] - z[iu[1]]
rng = np.random.default_rng(int(sys.argv[2]) if len(sys.argv) > 2 else 1)
t0 = time.time()
for s in [float(x) for x in sys.argv[1].split(',')]:
    best = (np.inf, None)
    for st in range(25):
        u0 = rng.uniform(-1, 1, d); u0[0] = 1.0
        cons = [dict(type='eq', fun=lambda u: np.array([u[0] - 1]), jac=lambda u: np.eye(d)[0][None, :])]
        if s > 0: cons.append(dict(type='ineq', fun=lambda u: gaps(u) ** 2 - s ** 2, jac=lambda u: 2 * gaps(u)[:, None] * Jz))
        sc = 1 / max(F(u0), 1e-300)
        r = minimize(lambda u: F(u) * sc, u0, jac=lambda u: G(u) * sc, constraints=cons, bounds=[(-1, 1)] * d, method='SLSQP', options=dict(maxiter=500, ftol=1e-18))
        u = r.x
        if np.abs(gaps(u)).min() >= s - 1e-7 and F(u) ** .5 < best[0]: best = (F(u) ** .5, u.copy())
    u = best[1]
    print('s>=%.3f: min |Phi|_2/(|c|_1 d^8) = %.3e  sep=%.4f sorted=%s (%.0fs)' % (s, best[0], np.abs(gaps(u)).min(), np.round(np.sort(np.concatenate([[0], u])), 4).tolist(), time.time() - t0), flush=True)
    np.save('sdc11_m8_land_s%.3f.npy' % s, u)
