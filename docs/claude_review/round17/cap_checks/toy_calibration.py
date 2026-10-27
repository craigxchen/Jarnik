# EVIDENCE ONLY (floating point). Toy sdc_q boxes A(q,L) = {k in {-1,0,1}^q : |k(P)| <= L, all 2q patterns}.
# For degree m: Psi_A(u) = min_P E_{mu_A}(l_u^m - P)^2 (L2 top component), Psi_T(u) exact cube formula
#   Psi_T(u) = (m!)^2 [t^m] prod_x (1 + (2/3) u_x^2 t + (1/18) u_x^4 t^2),
# restriction ratio rho_A = mu(A) Psi_A / Psi_T. Compare: Fourier circles u_theta (several theta, a = 1, 2),
# coordinate-permuted circles (same value multiset, structure destroyed), uniform random u. All u normalised |u|_inf = 1.
import sys, itertools, numpy as np, scipy.linalg as sla, math, time
from sdclib import patterns
def box(q, L):
    P, Pp = patterns(q)
    M = np.zeros((2 * q, q), dtype=np.int8)
    for i, S in enumerate(list(P) + list(Pp)):
        for x in S: M[i, x] = 1
    pts = []
    # enumerate cube in chunks
    allk = np.array(list(itertools.product((-1, 0, 1), repeat=q)), dtype=np.int8)
    s = allk.astype(np.int32) @ M.T.astype(np.int32)
    ok = np.all(np.abs(s) <= L, axis=1)
    return allk[ok]
def monos(q, mmax):
    out = []
    def rec(x, cur, deg):
        if x == q:
            out.append(tuple(cur)); return
        for a in range(3):
            if deg + a <= mmax:
                cur.append(a); rec(x + 1, cur, deg + a); cur.pop()
    rec(0, [], 0)
    return out
def psiT(u, m):
    poly = np.array([1.0])
    for ux in u:
        poly = np.convolve(poly, np.array([1.0, (2 / 3) * ux**2, ux**4 / 18]))[: m + 1]
    return math.factorial(m)**2 * poly[m] if len(poly) > m else 0.0
def run(q, L, ms, rng):
    t0 = time.time()
    A = box(q, L).astype(np.float64); U = len(A)
    print(f"q={q} L={L} |A|={U} mu(A)={U/3**q:.4f}  ({time.time()-t0:.1f}s)", flush=True)
    res = []
    for m in ms:
        mon = [a for a in monos(q, m - 1)]
        X = np.ones((U, len(mon)))
        for c, a in enumerate(mon):
            for x, e in enumerate(a):
                if e: X[:, c] *= A[:, x]**e
        Qm, R, piv = sla.qr(X, mode='economic', pivoting=True)
        dg = np.abs(np.diag(R)); rank = int(np.sum(dg > dg[0] * 1e-10))
        Qm = Qm[:, :rank]
        def psiA(u):
            f = (A @ u)**m
            r = f - Qm @ (Qm.T @ f)
            return float(r @ r) / U
        dirs = []
        for a in (1, 2):
            for th in (0.0, 0.37, 1.1, 2.0):
                u = np.cos(th + 2 * np.pi * a * np.arange(q) / q); dirs.append((f"circle a={a} th={th}", u))
        base = np.cos(0.37 + 2 * np.pi * np.arange(q) / q)
        for t in range(4):
            dirs.append((f"perm-circle #{t}", rng.permutation(base)))
        for t in range(4):
            dirs.append((f"uniform #{t}", rng.uniform(-1, 1, q)))
        print(f"  m={m}: #monomials={len(mon)} rank Poly_<m|_A={rank} dim V(A,m)={U-rank}", flush=True)
        for name, u in dirs:
            u = u / np.max(np.abs(u))
            pa, pt = psiA(u), psiT(u, m)
            rho = (U / 3**q) * pa / pt if pt > 0 else float('nan')
            line = f"    {name:24s} Psi_A={pa:.6e} Psi_T={pt:.6e} rho_A={rho:.5f}"
            print(line, flush=True); res.append((q, L, m, name, pa, pt, rho))
    return res
if __name__ == "__main__":
    rng = np.random.default_rng(17)
    cfg = [(7, 1, [2, 3, 4]), (7, 2, [3, 4, 5, 6, 7, 8]), (11, 1, [2, 3]), (11, 2, [3, 4, 5, 6, 7])]
    for q, L, ms in cfg:
        run(q, L, ms, rng)
