"""Numerical ZE test on a faithful model (evidence; exact verification of the rational subspaces found).
Input: json from had8.py (A, basis of V(A,m)).
(a) tangent cone: minimise tau(u) = |Q^T l_u^m|^2/|l_u^m|^2 over the unit sphere (Q = orthonormal basis of V);
(b) base locus near e: minimise G(u) = |Q^T g_u|^2/|u|^{2m}, g_u(k) = e^{i<k,u>} - sum_{t<m} (i<k,u>)^t/t!,
    on spheres |u| = eps;
(c) every zero is analysed by LLL integer-relation detection (relations v with |<v,u>| < tol, |v|_inf <= 50);
    degenerate iff at least one relation.
"""
import sys, json, math, time
import numpy as np
from scipy.optimize import minimize
import flint

def load(fn):
    D = json.load(open(fn))
    A = np.array(D['A'], dtype=float)
    B = np.array(D['basis'], dtype=float).T          # |A| x dimV
    Q, _ = np.linalg.qr(B)
    return D, A, Q

def tau(u, A, Q, m):
    x = A @ u
    f = x ** m
    nf = f @ f
    r = Q.T @ f
    return (r @ r) / nf

def tail(x, m):
    """e^{ix} - sum_{t<m} (ix)^t/t!  accurately (series for |x|<4, direct otherwise)"""
    out = np.zeros_like(x, dtype=complex)
    small = np.abs(x) < 4
    xs = x[small]
    term = (1j * xs) ** m / math.factorial(m)
    acc = term.copy()
    for t in range(m + 1, m + 60):
        term = term * (1j * xs) / t
        acc += term
    out[small] = acc
    xb = x[~small]
    s = np.zeros_like(xb, dtype=complex)
    for t in range(m):
        s += (1j * xb) ** t / math.factorial(t)
    out[~small] = np.exp(1j * xb) - s
    return out

def G(u, A, Q, m):
    g = tail(A @ u, m)
    r = Q.T @ g
    return float(np.real(np.vdot(r, r))) / (np.linalg.norm(u) ** (2 * m)) * math.factorial(m) ** 2

def relations(u, tol=1e-9, bound=60):
    """integer relations of u (LLL on [I | K u]); returns list of short v with |<v,u>| < tol*|v|."""
    d = len(u)
    K = 10 ** 12
    rows = []
    for i in range(d):
        rows.append([1 if j == i else 0 for j in range(d)] + [int(round(K * u[i]))])
    R = flint.fmpz_mat(rows).lll()
    rels = []
    for i in range(R.nrows()):
        v = [int(R[i, j]) for j in range(d)]
        if max(abs(x) for x in v) <= bound and abs(np.dot(v, u)) < tol * max(1, np.linalg.norm(v)):
            rels.append(v)
    return rels

def sphere_min(fun, d, starts, rng, scale=1.0):
    res = []
    for s in range(starts):
        u0 = rng.normal(size=d); u0 /= np.linalg.norm(u0)
        f = lambda z: fun(scale * z / np.linalg.norm(z))
        r = minimize(f, u0, method='BFGS', options=dict(gtol=1e-14, maxiter=2000))
        z = r.x / np.linalg.norm(r.x)
        res.append((fun(scale * z), z))
    return res

if __name__ == "__main__":
    fn = sys.argv[1]; starts = int(sys.argv[2]); mode = sys.argv[3]
    D, A, Q = load(fn)
    m = D['m']; d = A.shape[1]
    rng = np.random.default_rng(int(sys.argv[4]) if len(sys.argv) > 4 else 0)
    out = []
    def pr(*a):
        s = " ".join(str(x) for x in a); print(s, flush=True); out.append(s)
    pr("model", fn, "|A|=%d dimV=%d m=%d d=%d mode=%s starts=%d" % (A.shape[0], Q.shape[1], m, d, mode, starts))
    if mode == "tc":
        res = sphere_min(lambda u: tau(u, A, Q, m), d, starts, rng)
        thr = 1e-20
    else:
        eps = float(mode)
        res = sphere_min(lambda u: G(u, A, Q, m), d, starts, rng, scale=eps)
        thr = 1e-16
    vals = sorted(v for v, z in res)
    pr("minimum values: smallest 5", ["%.2e" % v for v in vals[:5]], "largest", "%.2e" % vals[-1])
    zeros = [(v, z) for v, z in res if v < thr]
    pr("zeros (value < %.0e): %d of %d" % (thr, len(zeros), len(res)))
    nondeg = 0
    ranks = {}
    for v, z in zeros:
        rel = relations(z)
        rk = np.linalg.matrix_rank(np.array(rel, dtype=float)) if rel else 0
        ranks[rk] = ranks.get(rk, 0) + 1
        if rk == 0:
            nondeg += 1
            pr("  NONDEGENERATE candidate:", "%.3e" % v, np.round(z, 6).tolist())
    pr("relation-lattice rank histogram of zeros (rank r => zero lies in a rational subspace of dim d-r):", ranks)
    pr("nondegenerate zeros:", nondeg)
    nonzeros = sorted([(v, z) for v, z in res if v >= thr], key=lambda t: t[0])
    if nonzeros:
        pr("smallest non-zero local minimum: %.3e" % nonzeros[0][0])
    open(fn.replace("_basis.json", "_zetest_%s_out.txt" % mode), "w").write("\n".join(out) + "\n")
