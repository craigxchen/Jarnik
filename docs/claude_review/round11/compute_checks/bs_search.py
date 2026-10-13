"""Base locus of a T = 1 model on the real torus (evidence, with exact coefficients and 50-digit refinement).
For k_x in {-1,0,1}: e^{i u_x k_x} = 1 + a_x k_x + b_x k_x^2, a_x = i sin u_x, b_x = cos u_x - 1.  Hence
f_u = sum_{alpha in {0,1,2}^d} c_alpha(u) k^alpha, c_alpha = prod_{alpha_x=1} a_x prod_{alpha_x=2} b_x, and
f_u in Poly_{<m}|_A  iff  sum_{|alpha|>=m} n_{j,alpha} c_alpha(u) = 0 for every basis vector c_j of V(A,m),
n_{j,alpha} = sum_k c_j(k) k^alpha (exact integers; the terms |alpha| < m drop out since c_j is in V(A,m)).
Trivial zeros (u with <= 3 nonzero coordinates when m = 8, more generally c(u) = 0) are removed by the
normalisation Psi(u) = |N c(u)|^2/|c(u)|^2.  We minimise Psi on spheres |u| = eps and on the whole torus,
refine zeros in mpmath, and classify them by support and integer relations."""
import sys, json, itertools, math
import numpy as np
from scipy.optimize import minimize
import mpmath as mp
import flint

fn = sys.argv[1]; starts = int(sys.argv[2]); eps = float(sys.argv[3]); seed = int(sys.argv[4])
D = json.load(open(fn))
A = D['A']; basis = D['basis']; m = D['m']; d = len(A[0])
alphas = [a for a in itertools.product((0, 1, 2), repeat=d) if sum(a) >= m]
Aarr = np.array(A, dtype=np.int64)
Bm = np.array(basis, dtype=object)
# exact n_{j,alpha}
Kpow = {}
cols = []
for a in alphas:
    v = np.ones(len(A), dtype=object)
    for x in range(d):
        if a[x]:
            v = v * (Aarr[:, x].astype(object) ** a[x])
    cols.append(v)
Nex = [[int(sum(c * v for c, v in zip(cj, col))) for col in cols] for cj in basis]
# sanity: low-degree moments vanish exactly (spot check)
low = [a for a in itertools.product((0, 1, 2), repeat=d) if sum(a) < m]
for a in low[:50]:
    v = np.ones(len(A), dtype=object)
    for x in range(d):
        if a[x]:
            v = v * (Aarr[:, x].astype(object) ** a[x])
    assert all(int(sum(c * w for c, w in zip(cj, v))) == 0 for cj in basis)
N = np.array(Nex, dtype=float)
Qn, _ = np.linalg.qr(N.T)
Nn = Qn.T
E = np.array(alphas)
ones = (E == 1); twos = (E == 2)

def cvec(u):
    a = 1j * np.sin(u); b = np.cos(u) - 1
    la = np.where(ones, a[None, :], 1.0); lb = np.where(twos, b[None, :], 1.0)
    return np.prod(la, axis=1) * np.prod(lb, axis=1)

def Psi(u):
    c = cvec(u)
    n2 = np.real(np.vdot(c, c))
    if n2 < 1e-300:
        return 1.0
    w = Nn @ c
    return float(np.real(np.vdot(w, w)) / n2)

out = []
def pr(*a):
    s = " ".join(str(x) for x in a); print(s, flush=True); out.append(s)
pr("model %s: d=%d m=%d dimV=%d, monomials |alpha|>=m: %d, eps=%s" % (fn, d, m, len(basis), len(alphas), eps))
rng = np.random.default_rng(seed)
zeros, vals = [], []
for s in range(starts):
    u0 = rng.normal(size=d); u0 /= np.linalg.norm(u0)
    if eps > 0:
        f = lambda z: Psi(eps * z / np.linalg.norm(z))
    else:
        u0 = rng.uniform(-math.pi, math.pi, size=d)
        f = Psi
    r = minimize(f, u0, method='BFGS', options=dict(gtol=1e-14, maxiter=4000))
    z = eps * r.x / np.linalg.norm(r.x) if eps > 0 else r.x
    v = Psi(z)
    vals.append(v)
    if v < 1e-16:
        zeros.append(z)
vals.sort()
pr("Psi local minima: smallest", ["%.1e" % v for v in vals[:4]], "median %.2e" % vals[len(vals) // 2], "#zeros(<1e-16): %d / %d" % (len(zeros), starts))
# refinement + classification
mp.mp.dps = 50
Nmp = [[mp.mpf(x) for x in row] for row in Nex]
def eqs(u):
    a = [mp.mpc(0, 1) * mp.sin(t) for t in u]; b = [mp.cos(t) - 1 for t in u]
    c = []
    for al in alphas:
        p = mp.mpc(1)
        for x in range(d):
            if al[x] == 1: p *= a[x]
            elif al[x] == 2: p *= b[x]
        c.append(p)
    nc = mp.sqrt(mp.fsum(abs(x) ** 2 for x in c))
    res = []
    for row in Nmp:
        s_ = mp.fsum(r * x for r, x in zip(row, c)) / nc
        res += [mp.re(s_), mp.im(s_)]
    return res
def refine(z):
    u = [mp.mpf(float(x)) for x in z]
    for it in range(30):
        res = eqs(u)
        if eps > 0:
            res.append(mp.fsum(t * t for t in u) - mp.mpf(eps) ** 2)
        h = mp.mpf(10) ** -25
        cols_ = []
        for t in range(d):
            up = list(u); up[t] += h
            rp = eqs(up)
            if eps > 0:
                rp.append(mp.fsum(x * x for x in up) - mp.mpf(eps) ** 2)
            cols_.append([(p - q) / h for p, q in zip(rp, res)])
        J = mp.matrix(len(res), d)
        for i in range(len(res)):
            for t in range(d):
                J[i, t] = cols_[t][i]
        rv = mp.matrix(res)
        step = mp.lu_solve(J.T * J + mp.eye(d) * mp.mpf(10) ** -45, J.T * rv)
        u = [u[t] - step[t] for t in range(d)]
        if mp.norm(rv) < mp.mpf(10) ** -42:
            break
    return u, mp.norm(mp.matrix(eqs(u)))
def relations(u, digits=35):
    K = mp.mpf(10) ** digits
    rows_ = [[1 if j == i else 0 for j in range(d)] + [int(mp.nint(K * u[i]))] for i in range(d)]
    R = flint.fmpz_mat(rows_).lll()
    rel = []
    for i in range(R.nrows()):
        v = [int(R[i, j]) for j in range(d)]
        if any(v) and max(abs(x) for x in v) <= 1000 and abs(mp.fsum(v[t] * u[t] for t in range(d))) < mp.mpf(10) ** -30:
            rel.append(v)
    return rel
classes = {}
nondeg = []
for z in zeros:
    u, res = refine(z)
    supp = tuple(t for t in range(d) if abs(u[t]) > mp.mpf(10) ** -20)
    rel = relations(u)
    rk = int(np.linalg.matrix_rank(np.array(rel, dtype=float))) if rel else 0
    key = (len(supp), rk, "res<1e-35" if res < mp.mpf(10) ** -35 else "res=%.0e" % float(res))
    classes[key] = classes.get(key, 0) + 1
    if rk == 0:
        nondeg.append((float(res), [float(x) for x in u]))
pr("refined zeros by (support size, relation rank, residual):", classes)
pr("NONDEGENERATE zeros:", len(nondeg))
for r_, u in nondeg[:5]:
    pr("   ", "%.2e" % r_, np.round(u, 10).tolist())
open(fn.replace("_basis.json", "_bs_eps%s_out.txt" % eps), "w").write("\n".join(out) + "\n")
