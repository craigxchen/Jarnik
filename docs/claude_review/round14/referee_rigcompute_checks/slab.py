# slab.py -- for a joint witness (per-unit demands m_(v,s), norm <= K), search for angles phi in R^r with
#   |sin(<v,phi> - s pi/4)| >= nu_s m_(v,s) / sqrt(Norm Z_v)   for all v in V_K, s in Z/4   (sharp.md (1.2)).
# A found phi is re-verified with mpmath at 50 digits (rigorous up to the stated margin).  Failure = evidence only.
import sys, math, random, pickle, json
import numpy as np
from scipy.optimize import minimize
import mpmath as mp
from reflib import *
from rverify import R
from repair import demands

def setup(P, K, labs):
    r = len(P); V = monos(r, K, half=True)
    dem = demands(P, K, labs, V)
    A = np.array(V, dtype=float); a = np.zeros((len(V), 4))
    for i, v in enumerate(V):
        N = norm_v(P, v)
        for s in range(4):
            a[i, s] = (1 if s % 2 == 0 else 2 ** -0.5) * dem.get((i, s), 1) / math.sqrt(N)
    return V, A, a, dem

def slack(phi, A, a):
    x = A @ phi
    vals = np.abs(np.sin(x[:, None] - np.arange(4)[None, :] * np.pi / 4)) - a
    return vals.min()

def search(A, a, r, tries=4000, seed=1):
    rng = np.random.default_rng(seed); best = (-9, None)
    for t in range(tries):
        phi = rng.uniform(0, np.pi, r)
        f = slack(phi, A, a)
        if f > best[0]: best = (f, phi)
        if f > 0: break
    # local polish
    res = minimize(lambda p: -slack(p, A, a), best[1], method="Nelder-Mead", options=dict(maxiter=4000, xatol=1e-12, fatol=1e-14))
    if -res.fun > best[0]: best = (-res.fun, res.x)
    return best

def verify_mp(phi, V, P, dem):
    mp.mp.dps = 50; worst = mp.mpf(10)
    for i, v in enumerate(V):
        x = mp.fsum(mp.mpf(t) * mp.mpf(p) for t, p in zip(v, phi)); N = norm_v(P, v)
        for s in range(4):
            lhs = abs(mp.sin(x - s * mp.pi / 4)); rhs = (1 if s % 2 == 0 else 1 / mp.sqrt(2)) * dem.get((i, s), 1) / mp.sqrt(N)
            worst = min(worst, lhs - rhs)
    return worst

if __name__ == "__main__":
    fn = sys.argv[1]
    W = R(open(fn, "rb")).load(); P, K = W["P"], W["K"]; r = len(P)
    labs = {q: [int(x) for x in W["labels"][q]["lifts"][-1]] for q in W["labels"]}
    V, A, a, dem = setup(P, K, labs)
    f, phi = search(A, a, r)
    w = verify_mp([float(x) for x in phi], V, P, dem)
    print(fn.split("/")[-1], "best slab slack (float)", f, "mp-verified worst slack", mp.nstr(w, 6), "phi", [float(x) for x in phi], flush=True)
    # genuine sanity: phi = arg pi_j with genuine demands on the same moduli
    pis = [gprime(p) for p in P]; ag = np.zeros_like(a)
    for i, v in enumerate(V):
        Dv = Dvals(Zv(pis, v)); N = norm_v(P, v)
        for t in range(4):
            g = 1
            for q in labs:
                e = val(Dv[t], q)
                if e: g *= q ** e
            ag[i, t] = (1 if t % 2 == 0 else 2 ** -0.5) * g / math.sqrt(N)
    phig = np.array([math.atan2(b, a_) for a_, b in pis])
    print("   genuine demand on the same moduli, phi = arg pi_j: slack", slack(phig, A, ag), flush=True)
