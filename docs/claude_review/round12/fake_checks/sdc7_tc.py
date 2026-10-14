"""sdc7 model (T=1, L=2, m=8, dim V=1): the single tangent form Phi = sum_k c_k <k,u>^8.
(1) exact c (recomputed), parity, |c|_1;
(2) max over the real zero set {Phi = 0} of the minimal pair separation  sep(u) = min_{x<y in rows} |u_x - u_y|
    (base row u_0 = 0), normalised |u|_inf = 1;  and the min of Gamma(u) = (2/N2) sum ln^+(1/|u_x-u_y|) on {Phi=0};
(3) compare with the same quantities over all u (equally spaced optimum) -> is the zero set 'comfortable' for fakes?
(4) sign structure: max Phi/(|c|_1 d^8) over {sep >= s} for a grid of s (where is Phi > 0?).
Numerics are evidence; the certified witness is produced by sdc7_cert.py.
"""
import json, time, sys
import numpy as np
from scipy.optimize import minimize
from flib import *
from model import V_basis

M, pats = sdc_patterns(7)
A = box(M, pats, 2)
d = M - 1
t = time.time()
B = V_basis(A, d, 8)
assert len(B) == 1
c = np.array(B[0], dtype=np.int64)
if c.sum() < 0 or (c[np.nonzero(c)[0][0]] < 0):
    pass
Aa = np.array(A, dtype=np.int64)
idx = {tuple(k): i for i, k in enumerate(A)}
even = all(c[idx[tuple(-np.array(k))]] == c[i] for i, k in enumerate(A))
c1 = int(np.abs(c).sum())
print('|A|=%d dim V=1, |c|_inf=%d |c|_1=%d support=%d even=%s (%.1fs)' % (len(A), np.abs(c).max(), c1, np.count_nonzero(c), even, time.time() - t))
json.dump(dict(A=A, c=[int(x) for x in c]), open('sdc7_c.json', 'w'))
cf = c.astype(float)
Af = Aa.astype(float)


def Phi(u):
    return float(cf @ (Af @ u) ** 8)


def gradPhi(u):
    y = Af @ u
    return 8 * (cf * y ** 7) @ Af


def kappa(u):
    return abs(Phi(u)) / (c1 * d ** 8)


def pairs(u):
    z = np.concatenate([[0.0], u])
    return np.abs(z[:, None] - z[None, :])[np.triu_indices(len(z), 1)]


def sep(u):
    return pairs(u).min()


def Gamma(u):
    pz = pairs(u)
    N2 = len(pz)
    return 2.0 / N2 * np.sum(np.log(np.maximum(1.0, 1.0 / np.maximum(pz, 1e-300))))


rng = np.random.default_rng(1)
# (2) maximise s subject to Phi(u)=0 (scaled), -1<=u<=1, pair seps >= s.  variables (u, s)
iu = np.triu_indices(d + 1, 1)


def cons_pairs(x):
    u, s = x[:-1], x[-1]
    z = np.concatenate([[0.0], u])
    D = z[iu[0]] - z[iu[1]]
    return np.concatenate([D ** 2 - s ** 2])


best = None
for trial in range(300):
    u0 = rng.uniform(-1, 1, d)
    x0 = np.concatenate([u0, [0.01]])
    cons = [dict(type='eq', fun=lambda x: np.array([Phi(x[:-1]) / (c1 * 1.0)]),
                 jac=lambda x: np.concatenate([gradPhi(x[:-1]) / c1, [0.0]])[None, :]),
            dict(type='ineq', fun=cons_pairs)]
    res = minimize(lambda x: -x[-1], x0, jac=lambda x: np.concatenate([np.zeros(d), [-1.0]]),
                   constraints=cons, bounds=[(-1, 1)] * d + [(0, 1)], method='SLSQP', options=dict(maxiter=500, ftol=1e-14))
    u = res.x[:-1]
    if abs(Phi(u)) / c1 < 1e-9 and np.abs(u).max() > 0.999:
        s = sep(u)
        if best is None or s > best[0]:
            best = (s, u.copy())
print('max sep on {Phi=0} (|u|_inf=1): %.6f at u=%s, Gamma=%.4f' % (best[0], np.round(best[1], 6).tolist(), Gamma(best[1])))
np.save('sdc7_bestsep.npy', best[1])
# (3) unconstrained: equally spaced points have sep 2/7
print('max sep over all u with |u|_inf=1: 2/7 = %.6f (equally spaced); its Phi sign: %+.3e' % (
    2 / 7, Phi(np.linspace(-1, 1, 8)[1:] - np.linspace(-1, 1, 8)[0] - 1) / c1))
# (4) on {sep >= s}: maximise and minimise Phi/(|c|_1 d^8)
for s in [0.0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.28]:
    vals = []
    for sign in (+1, -1):
        bestv = None
        for trial in range(60):
            u0 = rng.uniform(-1, 1, d)
            cons = [dict(type='ineq', fun=lambda u, s=s: (lambda z: (z[iu[0]] - z[iu[1]]) ** 2 - s ** 2)(np.concatenate([[0.0], u])))]
            res = minimize(lambda u: -sign * Phi(u) / (c1 * d ** 8), u0, jac=lambda u: -sign * gradPhi(u) / (c1 * d ** 8),
                           constraints=cons, bounds=[(-1, 1)] * d, method='SLSQP', options=dict(maxiter=300))
            if res.success or True:
                if sep(res.x) >= s - 1e-6:
                    v = Phi(res.x) / (c1 * d ** 8)
                    if bestv is None or sign * v > sign * bestv:
                        bestv = v
        vals.append(bestv)
    print('sep >= %.2f: max kappa-signed Phi/(|c|_1 d^8) = %+.3e, min = %+.3e' % (s, vals[0], vals[1]))
