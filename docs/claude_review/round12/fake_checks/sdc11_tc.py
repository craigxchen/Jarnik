"""sdc11, L=2, m=8: common real zeros of the 5 tangent octics (evidence; certified afterwards by sdc11_cert.py).
Maximise the separation s over {Phi_i(u) = 0 (i<=5), u_0 = 1, |u|_inf <= 1, all pair gaps >= s} (base row = 0).
The forms are invariant under the cyclic shift of the 11 coordinates (c is shift-invariant), so u_0 = 1 is WLOG.
"""
import json, time
import numpy as np
from scipy.optimize import minimize

D = json.load(open('sdc11_m8_basis.json'))
A = np.array(D['A'], float)
C = np.array(D['basis'], float)
C = C / np.abs(C).sum(axis=1, keepdims=True)  # |c|_1 = 1
d = 11; m = 8


def Phis(u):
    return C @ (A @ u) ** m


def J(u):
    y = A @ u
    return m * (C * y ** (m - 1)) @ A


iu = np.triu_indices(d + 1, 1)


def gaps(u):
    z = np.concatenate([[0.0], u])
    return z[iu[0]] - z[iu[1]]


Jz = np.zeros((len(iu[0]), d + 1))
Jz[np.arange(len(iu[0])), iu[0]] = 1
Jz[np.arange(len(iu[0])), iu[1]] -= 1
Jz = Jz[:, 1:]
rng = np.random.default_rng(3)
best = []
t0 = time.time()
NST = 60
for trial in range(NST):
    u0 = rng.uniform(-1, 1, d); u0[0] = 1.0
    x0 = np.concatenate([u0, [0.005]])
    cons = [dict(type='eq', fun=lambda x: Phis(x[:-1]) * d ** m, jac=lambda x: np.hstack([J(x[:-1]) * d ** m, np.zeros((5, 1))])),
            dict(type='eq', fun=lambda x: np.array([x[0] - 1.0]), jac=lambda x: np.eye(d + 1)[0][None, :]),
            dict(type='ineq', fun=lambda x: gaps(x[:-1]) ** 2 - x[-1] ** 2,
                 jac=lambda x: np.hstack([2 * gaps(x[:-1])[:, None] * Jz, -2 * x[-1] * np.ones((len(iu[0]), 1))]))]
    res = minimize(lambda x: -x[-1], x0, jac=lambda x: np.concatenate([np.zeros(d), [-1.0]]), constraints=cons,
                   bounds=[(-1, 1)] * d + [(0, 1)], method='SLSQP', options=dict(maxiter=1000, ftol=1e-15))
    u = res.x[:-1]
    resid = np.abs(Phis(u)).max() * d ** m
    s = np.abs(gaps(u)).min()
    if resid < 1e-9 and abs(u[0] - 1) < 1e-9:
        sv = np.linalg.svd(J(u), compute_uv=False)
        best.append((s, u.copy(), sv[-1] / sv[0]))
        print('trial %d: zero with sep=%.5f cond=%.2e (%.0fs)' % (trial, s, sv[-1]/sv[0], time.time()-t0), flush=True)
        np.save('sdc11_zeros.npy', np.array([b[1] for b in sorted(best, key=lambda t: -t[0])]))
best.sort(key=lambda t: -t[0])
print('%d converged zeros out of %d starts (%.0fs)' % (len(best), NST, time.time() - t0))
for s, u, cond in best[:8]:
    print('sep=%.5f  Jacobian sigma_min/sigma_max=%.2e  sorted=%s' % (s, cond, np.round(np.sort(np.concatenate([[0], u])), 4).tolist()))
print('equally spaced optimum: 2/11 = %.5f' % (2 / 11))
np.save('sdc11_zeros.npy', np.array([b[1] for b in best]))
