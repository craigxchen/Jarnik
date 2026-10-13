# EVIDENCE. Single-J tangent cones for supercritical small models: are all real TC points rational (degenerate)?
import sys, numpy as np, random, math
from scipy.optimize import minimize
from engine import *
seed = int(sys.argv[1]) if len(sys.argv) > 1 else 1
rng = np.random.default_rng(seed); random.seed(seed)
def rand_profile(n, r):
    cols = []
    while len(cols) < r:
        a = tuple(random.randint(0, 1) for _ in range(n))
        if 0 < sum(a) < n: cols.append(a)
    return cols
for (n, r, T, w) in [(4, 4, 4, 6), (4, 8, 4, 6), (5, 5, 3, 6), (5, 10, 3, 6)]:
    cols = rand_profile(n, r); J = list(range(n)); base = n-1
    omega = [1.0/r]*r
    m, U, dV, lp, mh = choose_m(J, base, cols, T, w, omega, 1)
    print(f"\nn={n} r={r} T={T} w={w}: |A|={U}, LP sum omega h = {lp:.2f}, max h = {mh}, m={m}, dim V(A,m)={dV}", flush=True)
    if dV < 2*(n-1): print("  dim V too small; skip"); continue
    S = System(J, base, cols, T, w, m)
    print(f"  rank Poly={S.rk}, singular gap {S.gap[0]:.1e} | {S.gap[1]:.1e}", flush=True)
    def obj(u):
        u = np.append(u, 0.0); u = u - u.mean(); u /= np.linalg.norm(u); return S.tau(u)
    sols = []
    for s in range(200):
        r_ = minimize(obj, rng.normal(size=n-1), method='Nelder-Mead', options={'xatol':1e-13, 'fatol':1e-26, 'maxiter':20000})
        u = np.append(r_.x, 0.0); u = u - u.mean(); u /= np.linalg.norm(u)
        if u[np.argmax(np.abs(u))] < 0: u = -u
        sols.append((r_.fun, u))
    sols.sort(key=lambda t: t[0])
    zs = []
    for f, u in sols:
        if f < 1e-14 and not any(np.linalg.norm(u-z) < 1e-5 for z in zs): zs.append(u)
    print(f"  min tau = {sols[0][0]:.2e}; #distinct TC zeros (tau<1e-14): {len(zs)}; 10th smallest tau {sols[min(9,len(sols)-1)][0]:.2e}")
    for z in zs:
        # rationality: express in basis; find integer relation pattern
        rel = int_relation(z, vmax=4, tol=1e-6)
        rat = z/np.abs(z[np.abs(z) > 1e-6]).min()
        print("   TC zero", np.round(rat, 5), " smallest relation", rel)
