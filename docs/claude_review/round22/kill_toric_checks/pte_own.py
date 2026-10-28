#!/usr/bin/env python3
"""Independent numerical search (evidence only) for toric.md §7: balanced V-null rational curves of
minimal degree (Gaussian linear blocks), k points.  Unknowns r_T (T nonempty subset of [k-1]),
gauge r_{[k-1]} = i.  Equations: Im of coeff of s^1..s^{deg-1} of P_j = prod_{T∋j}(s - r_T), j in [k-1].
Own code; scipy least_squares (trf) with random starts of several scales.  usage: pte_own.py k nstarts seed"""
import sys, itertools
import numpy as np
from scipy.optimize import least_squares
k = int(sys.argv[1]); nst = int(sys.argv[2]); seed = int(sys.argv[3])
rng = np.random.default_rng(seed)
n = k-1
blocks = [T for r in range(1, n+1) for T in itertools.combinations(range(n), r)]
top = tuple(range(n))
free = [T for T in blocks if T != top]
idx = {T: i for i, T in enumerate(free)}

def roots_of(z):
    r = {top: 1j}
    for T in free:
        r[T] = z[2*idx[T]] + 1j*z[2*idx[T]+1]
    return r

def polys(r):
    out = []
    for j in range(n):
        R = [r[T] for T in blocks if j in T]
        out.append(np.poly(R))   # leading coeff first
    return out

def resid(z):
    r = roots_of(z)
    res = []
    for p in polys(r):
        # p[0] = 1 (s^deg), p[-1] = constant; need Im of p[1..deg-1]
        res.extend(np.imag(p[1:-1]) / (1.0 + np.abs(p[1:-1])))
    return np.array(res)

def degenerate_info(z):
    r = roots_of(z)
    ps = polys(r)
    c = np.array([np.imag(p[-1]) for p in ps])
    scale = np.array([np.max(np.abs(p)) for p in ps])
    imr = min(abs(np.imag(v)) for v in r.values())
    allr = list(r.values()) + [np.conj(v) for v in r.values()]
    mind = min(abs(a-b) for a, b in itertools.combinations(allr, 2))
    cn = c/scale
    cdist = min([abs(cn[a]-cn[b]) for a, b in itertools.combinations(range(n), 2)] + [1.0])
    return cn, imr, mind, cdist

found = 0; conv = 0; stats = []
for st in range(nst):
    sc = rng.choice([0.5, 1.0, 2.0, 4.0])
    z0 = rng.normal(size=2*len(free))*sc
    try:
        sol = least_squares(resid, z0, method=('lm' if 2*len(free) <= n*(2**(n-1)-1) else 'trf'), xtol=1e-15, ftol=1e-15, gtol=1e-15, max_nfev=20000)
    except Exception as e:
        continue
    f = np.max(np.abs(sol.fun))
    if f < 1e-10:
        conv += 1
        cn, imr, mind, cdist = degenerate_info(sol.x)
        nondeg = (np.min(np.abs(cn)) > 1e-6) and imr > 1e-4 and mind > 1e-4 and cdist > 1e-6
        stats.append((np.min(np.abs(cn)), imr, mind))
        if nondeg:
            found += 1
            if found <= 3:
                print("NONDEGENERATE solution, start", st, "c/scale =", cn, "min|Im r| =", imr, "min root dist =", mind)
print(f"k={k} starts={nst} converged={conv} nondegenerate={found}")
if stats:
    a = np.array(stats)
    print("converged: min|c|/scale quantiles", np.quantile(a[:, 0], [0, .5, 1]), "min|Im r| quantiles", np.quantile(a[:, 1], [0, .5, 1]),
          "min root distance quantiles", np.quantile(a[:, 2], [0, .5, 1]))
