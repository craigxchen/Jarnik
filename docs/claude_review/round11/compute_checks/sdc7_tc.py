"""Real zeros of the tangent-cone form Phi (sdc7 model; Phi <= 0 on all tested integer points).
Minimise -Phi(u)/|u|^8 on the sphere from random starts (exact integer form coefficients, float evaluation),
then at each zero: Hessian rank of -Phi (local codimension of the zero set if it is a smooth Bott-type minimum),
and integer relations (LLL).  Evidence."""
import json, math, itertools
import numpy as np
from scipy.optimize import minimize
import flint
from core import monomials
D = json.load(open("sdc7_T1_L2_m8_basis.json")); A = np.array(D['A'], float); c = np.array(D['basis'][0], float); d = 7
mons = monomials(d, 8)
E = np.array(mons)
# exact integer coefficients of Phi/8! (python ints), then float
Ai = [tuple(int(x) for x in k) for k in D['A']]; ci = D['basis'][0]
coef = []
for a in mons:
    mult = math.factorial(8)
    for ai in a: mult //= math.factorial(ai)
    s = 0
    for k, ck in zip(Ai, ci):
        if ck:
            p = 1
            for t in range(d):
                if a[t]: p *= k[t] ** a[t]
            s += ck * p
    coef.append(mult * s // math.factorial(8))
    assert mult * s % math.factorial(8) == 0
print("Phi/8! coefficients: nonzero", sum(1 for x in coef if x), "values", sorted(set(coef)))
coef = np.array(coef, dtype=float)
def Phi(u):
    return coef @ np.prod(u[None, :] ** E, axis=1)
def f(u):
    u = u / np.linalg.norm(u)
    return -Phi(u)
rng = np.random.default_rng(1)
res = []
for s in range(300):
    u0 = rng.normal(size=d)
    r = minimize(f, u0, method='BFGS', options=dict(gtol=1e-12, maxiter=5000))
    u = r.x / np.linalg.norm(r.x)
    res.append((f(u), u))
vals = np.array([v for v, u in res])
print("min of -Phi on S^6 from 300 starts: smallest %.2e, largest local min %.2e, #(<1e-9): %d" % (vals.min(), vals.max(), (vals < 1e-9).sum()))
def hess_rank(u):
    h = 1e-4
    H = np.zeros((d, d))
    g = lambda z: -Phi(z)
    for i in range(d):
        for j in range(d):
            ei = np.eye(d)[i] * h; ej = np.eye(d)[j] * h
            H[i, j] = (g(u + ei + ej) - g(u + ei - ej) - g(u - ei + ej) + g(u - ei - ej)) / (4 * h * h)
    sv = np.linalg.svd(H, compute_uv=False)
    return int(np.sum(sv > 1e-6 * max(1, sv[0]))), sv
def relations(u, tol=1e-7):
    K = 10 ** 10
    rows = [[1 if j == i else 0 for j in range(d)] + [int(round(K * u[i]))] for i in range(d)]
    R = flint.fmpz_mat(rows).lll()
    rel = [[int(R[i, j]) for j in range(d)] for i in range(R.nrows())]
    rel = [v for v in rel if any(v) and max(abs(x) for x in v) <= 20 and abs(np.dot(v, u)) < tol]
    return int(np.linalg.matrix_rank(np.array(rel, float))) if rel else 0
stats = {}
for v, u in res:
    if v < 1e-9:
        hr, sv = hess_rank(u)
        rk = relations(u)
        key = (hr, 7 - rk)
        stats[key] = stats.get(key, 0) + 1
print("zeros by (Hessian rank of -Phi, dim of smallest rational subspace containing the zero):", stats)

# exact positivity witness: maximise Phi on the sphere, round to integers (scale 10^4), evaluate exactly
from fractions import Fraction
best = min(res, key=lambda t: t[0])
u = best[1]
for s in range(50):
    u0 = rng.normal(size=d)
    r = minimize(f, u0, method='BFGS', options=dict(gtol=1e-12, maxiter=5000))
    if f(r.x) < f(u):
        u = r.x / np.linalg.norm(r.x)
ui = [int(round(10 ** 4 * x)) for x in u]
coef_i = [int(x) for x in coef]
val = sum(cf * math.prod(ui[t] ** int(E[j][t]) for t in range(d)) for j, cf in enumerate(coef_i) if cf)
print("EXACT: Phi/8! at integer point u+ =", ui, "equals", val, "(> 0)" if val > 0 else "(not > 0)")
um = [1, 0, 0, 0, 0, 0, 0]
for cand in itertools.product(range(-1, 2), repeat=d):
    vv = sum(cf * math.prod(cand[t] ** int(E[j][t]) for t in range(d)) for j, cf in enumerate(coef_i) if cf)
    if vv < 0:
        print("EXACT: Phi/8! at u- =", list(cand), "equals", vv); break
