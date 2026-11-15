"""Function-field circle problem over R(t) (k = R, k(i) = C), squarefree N with roots a_k, abar_k.
Monic point z_eps = prod_k (t - a_k^{eps_k}), a^+ = a, a^- = conj(a).  deg(z_eps - z_eps') <= s := floor(D/2 + c)
iff e_r agree for r = 1..m, m = D - 1 - s, iff Im p_r agree for r = 1..m (Re parts agree automatically; Newton).
Question: for M sign vectors, does a non-degenerate a in C^D exist?  Generic count: 2D - 2 real unknowns (real affine
maps t -> alpha t + beta preserve everything) against (M-1) m equations.
Numerical test (EVIDENCE ONLY): Levenberg-Marquardt from random starts, random designs with pairwise distance >= m+1
(necessary by the PTE/Liouville lemma).  Normalisation: sum Re a = 0, sum |a|^2 = D.
"""
import sys, itertools
import numpy as np
from scipy.optimize import least_squares
rng = np.random.default_rng(int(sys.argv[4]) if len(sys.argv) > 4 else 1)
D = int(sys.argv[1]); M = int(sys.argv[2]); c = 0
ndes = int(sys.argv[3])
s = int(np.floor(D / 2 + c)); m = D - 1 - s
print(f"D={D} M={M} m={m}: unknowns {2*D-2}, equations {(M-1)*m}", flush=True)

def amap(x): return x[:D] + 1j * (0.1 + x[D:] ** 2)   # Im a_k >= 0.1 (a_k labels the upper root of each pair)
def resid(x, E):
    a = amap(x)
    W = np.array([np.imag(a ** r) for r in range(1, m + 1)])  # m x D
    V = W @ E.T  # m x M
    r = (V[:, 1:] - V[:, :1]).ravel()
    return np.concatenate([r, [np.sum(x[:D]), np.sum(np.abs(a) ** 2) - D]])

def nondeg(x):
    a = amap(x)
    d1 = np.min(np.abs(a.imag))
    d2 = min(min(abs(a[i] - a[j]), abs(a[i] - np.conj(a[j]))) for i in range(D) for j in range(i + 1, D))
    return min(d1, d2)

allE = np.array(list(itertools.product([1, -1], repeat=D)))
found = 0; bestres = 1e9; allbest = []
for des in range(ndes):
    # random design with pairwise Hamming distance >= m+1
    for _ in range(1000):
        idx = rng.choice(len(allE), size=M, replace=False); E = allE[idx]
        H = (D - E @ E.T) / 2
        if np.all(H[np.triu_indices(M, 1)] >= m + 1): break
    else:
        continue
    for start in range(40):
        x0 = rng.normal(size=2 * D)
        sol = least_squares(resid, x0, args=(E,), method='trf', xtol=1e-15, ftol=1e-15, gtol=1e-15, max_nfev=3000)
        r = np.max(np.abs(sol.fun)); nd = nondeg(sol.x)
        if nd > 0.05: bestres = min(bestres, r)
        allbest.append((r, nd))
        if r < 1e-10 and nd > 0.05:
            found += 1
            if found <= 3:
                print(" solution: design", E.tolist(), "resid", r, "nondeg", round(nd, 3), "a=", np.round(amap(sol.x), 4).tolist(), flush=True)
            break
print(f"designs tried {ndes}; designs with a nondegenerate solution: {found}; best nondegenerate residual {bestres:.2e}")
allbest.sort(); print("lowest residuals (res, nondeg):", [(f"{a:.1e}", round(b,3)) for a,b in allbest[:5]])
