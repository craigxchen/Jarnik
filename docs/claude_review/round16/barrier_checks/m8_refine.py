"""High-precision refinement of numerical TC zeros of Syl(8) boxes (M8: L=2,m=8; M8': L=3,m=11).
Gauss-Newton in mpmath (60 digits) on the overdetermined system {Phi_c(u)=0, c in exact basis of V} + {u_x0 = 1}.
A candidate is reported as a GENUINE collision-free zero only if the residual drops below 1e-45 while the minimal
gap among (0,u_1..u_7) stays >= 1e-3.  Then PSLQ looks for integer relations with small coefficients.
EVIDENCE-level (floating point); exact basis.  Usage: python3 m8_refine.py L m nstarts
"""
import sys, time
import numpy as np
import mpmath as mp
import flint
from flib import sylvester_patterns, box, reduced_monomials
from scipy.optimize import least_squares

L = int(sys.argv[1]); m = int(sys.argv[2]); NST = int(sys.argv[3])
mp.mp.dps = 60
M, pats = sylvester_patterns(3)
A = box(M, pats, L)
d = M - 1
mons = [mm for deg in range(m) for mm in reduced_monomials(d, deg)]
E = flint.fmpz_mat([[int(np.prod([k[i] ** e[i] for i in range(d)])) for k in A] for e in mons])
X, nul = E.nullspace()
basis = [[int(X[i, j]) for i in range(X.nrows())] for j in range(nul)]
print(f"Syl(8) L={L} m={m} |A|={len(A)} dimV={nul}", flush=True)
Amat = np.array(A, dtype=float)
C = np.array(basis, dtype=float); C /= np.abs(C).sum(axis=1, keepdims=True)
# integer forms for mp evaluation: Phi_c(u) = sum_k c_k (k.u)^m ; precompute supports
supp = [[(A[t], basis[j][t]) for t in range(len(A)) if basis[j][t] != 0] for j in range(nul)]
norm1 = [sum(abs(v) for _, v in s) for s in supp]


def res_np(u):
    v = Amat @ u
    return C @ (v ** m) / d ** m


def jac_np(u):
    v = Amat @ u
    return (C * (m * v ** (m - 1))) @ Amat / d ** m


def mp_res_jac(u):
    R = []; J = []
    for j in range(nul):
        r = mp.mpf(0); g = [mp.mpf(0)] * d
        for k, cv in supp[j]:
            s = mp.fsum(k[i] * u[i] for i in range(d) if k[i])
            r += cv * s ** m
            if s != 0:
                ds = cv * m * s ** (m - 1)
                for i in range(d):
                    if k[i]:
                        g[i] += ds * k[i]
        R.append(r / norm1[j]); J.append([gi / norm1[j] for gi in g])
    return R, J


rng = np.random.default_rng(7)
genuine = []
stats = {"slid_to_stratum": 0, "no_conv": 0, "genuine": 0}
t0 = time.time()
for s in range(NST):
    u0 = rng.uniform(-1, 1, d)
    f = lambda y: np.concatenate([res_np(y) * 1e6, [np.dot(y, y) - 1.0]])
    Jf = lambda y: np.vstack([jac_np(y) * 1e6, 2 * y[None, :]])
    r = least_squares(f, u0 / np.linalg.norm(u0), jac=Jf, xtol=1e-15, ftol=1e-15, gtol=1e-15, max_nfev=3000)
    u = r.x / np.abs(r.x).max()
    if np.abs(res_np(u)).max() > 1e-12:
        stats["no_conv"] += 1; continue
    x0 = int(np.argmax(np.abs(u)))
    um = [mp.mpf(float(v)) for v in u]
    um = [v / um[x0] for v in um]
    ok = False
    for it in range(40):
        R, J = mp_res_jac(um)
        # drop column x0 (fixed = 1)
        Jr = mp.matrix([[J[a][b] for b in range(d) if b != x0] for a in range(nul)])
        Rv = mp.matrix(R)
        try:
            delta = mp.lu_solve(Jr.T * Jr, -(Jr.T * Rv))
        except ZeroDivisionError:
            break
        idx = 0
        for b in range(d):
            if b != x0:
                um[b] += delta[idx]; idx += 1
        rn = max(abs(v) for v in R)
        z = [mp.mpf(0)] + um
        gap = min(abs(z[a] - z[b]) for a in range(M) for b in range(a + 1, M))
        if gap < mp.mpf('1e-4'):
            break
        if rn < mp.mpf('1e-45'):
            ok = True; break
    z = [mp.mpf(0)] + um
    gap = min(abs(z[a] - z[b]) for a in range(M) for b in range(a + 1, M))
    if ok and gap >= mp.mpf('1e-3'):
        stats["genuine"] += 1
        genuine.append((float(gap), [float(v) for v in um]))
        rel = mp.pslq(z[1:], maxcoeff=20, maxsteps=10 ** 5)
        print(f"  genuine zero #{stats['genuine']}: mingap={float(gap):.4f} u={[round(float(v),6) for v in um]} pslq={rel}", flush=True)
    else:
        stats["slid_to_stratum"] += 1
print(stats, f"{time.time()-t0:.0f}s")
