"""Base locus of a T = 1 model on the real torus near e, as a real polynomial system (evidence).
For k in {-1,0,1}: e^{iuk} = 1 + a k + b k^2 with a = i sin u, b = cos u - 1.  With tau = tan(u/2):
a = 2 i tau/(1+tau^2), b = -2 tau^2/(1+tau^2).  Clearing prod_x (1+tau_x^2):
   E_j(tau) = sum_{|alpha|>=m} n_{j,alpha} (2i)^{#1(alpha)} (-2)^{#2(alpha)} tau^alpha prod_{x notin supp alpha} (1+tau_x^2)
with exact integers n_{j,alpha} = sum_k c_j(k) k^alpha.  E_j is real (even c_j) or i*real (odd c_j) up to the
power of i.  We solve min sum_j (E_j/scale_j)^2 with scale_j = sum_alpha |terms| (relative residual), on spheres
|tau| = eps, from random starts, by Levenberg-Marquardt; zeros are classified by support and integer relations
(LLL), and the local dimension of the zero set is estimated from the numerical rank of the Jacobian."""
import sys, json, itertools, math
import numpy as np
from scipy.optimize import least_squares
import flint

fn = sys.argv[1]; starts = int(sys.argv[2]); eps = float(sys.argv[3]); seed = int(sys.argv[4])
D = json.load(open(fn))
A = D['A']; basis = D['basis']; m = D['m']; d = len(A[0])
alphas = [a for a in itertools.product((0, 1, 2), repeat=d) if sum(a) >= m]
Aarr = np.array(A, dtype=object)
cols = []
for a in alphas:
    v = np.ones(len(A), dtype=object)
    for x in range(d):
        if a[x]:
            v = v * (Aarr[:, x] ** a[x])
    cols.append(v)
Bo = np.array(basis, dtype=object)
Nex = np.array([[int(np.dot(cj, col)) for col in cols] for cj in Bo], dtype=float)   # dimV x #alpha
E = np.array(alphas); n1 = (E == 1).sum(1); n2 = (E == 2).sum(1)
# real coefficient: (2i)^{n1} (-2)^{n2} = i^{n1} 2^{n1+n2} (-1)^{n2}; equation j is real iff all its alphas have n1 even
coef = (2.0 ** (n1 + n2)) * ((-1.0) ** n2) * np.where(n1 % 4 < 2, 1.0, -1.0)   # i^{n1} = (+-1) * i^{n1 mod 2}
par = n1 % 2
supp = (E > 0)
def terms(t):
    t2 = 1 + t * t
    mon = np.prod(np.where(E > 0, t[None, :] ** E, 1.0), axis=1)
    rest = np.prod(np.where(supp, 1.0, t2[None, :]), axis=1)
    return coef * mon * rest           # times i^{par}
def resid(t):
    T_ = terms(t)
    out = []
    for p in (0, 1):
        mask = par == p
        num = Nex[:, mask] @ T_[mask]
        den = np.abs(Nex[:, mask]) @ np.abs(T_[mask]) + 1e-300
        out.append(num / den)
    return np.concatenate(out)
def relations(u, tol):
    K = 10 ** 11
    rows = [[1 if j == i else 0 for j in range(d)] + [int(round(K * u[i]))] for i in range(d)]
    R = flint.fmpz_mat(rows).lll()
    rel = []
    for i in range(R.nrows()):
        v = [int(R[i, j]) for j in range(d)]
        if any(v) and max(abs(x) for x in v) <= 30 and abs(np.dot(v, u)) < tol:
            rel.append(v)
    return rel
out = []
def pr(*a):
    s = " ".join(str(x) for x in a); print(s, flush=True); out.append(s)
pr("model %s eps=%g starts=%d (|alpha|>=m monomials: %d)" % (fn, eps, starts, len(alphas)))
rng = np.random.default_rng(seed)
classes = {}
cand = []
vals = []
for s in range(starts):
    z0 = rng.normal(size=d); z0 *= eps / np.linalg.norm(z0)
    f = lambda z: np.concatenate([resid(z), [(z @ z - eps * eps) / eps ** 2]])
    r = least_squares(f, z0, method='lm', xtol=1e-15, ftol=1e-15, gtol=1e-15, max_nfev=20000)
    z = r.x
    val = np.linalg.norm(resid(z))
    vals.append(val)
    if val > 1e-9:
        classes[("nonzero",)] = classes.get(("nonzero",), 0) + 1
        continue
    u = 2 * np.arctan(z)
    sup = tuple(i for i in range(d) if abs(u[i]) > 1e-7 * eps)
    rel = relations(u, 1e-6 * eps)
    rk = int(np.linalg.matrix_rank(np.array(rel, dtype=float))) if rel else 0
    # numerical local dimension: corank of Jacobian of resid at z
    h = 1e-7 * eps
    J = np.array([(resid(z + h * e) - resid(z - h * e)) / (2 * h) for e in np.eye(d)]).T
    sv = np.linalg.svd(J, compute_uv=False)
    loc = int(np.sum(sv < 1e-6 * sv[0]))
    key = ("zero", len(sup), d - rk, loc)
    classes[key] = classes.get(key, 0) + 1
    if rk == 0:
        cand.append(u)
vals.sort()
pr("residual minima: smallest %s ; largest %.2e" % (["%.1e" % v for v in vals[:3]], vals[-1]))
pr("classes (zero, support size, dim of smallest rational subspace containing it, numerical corank of Jacobian):")
for k in sorted(classes, key=str):
    pr("   ", k, classes[k])
pr("zeros with NO integer relation (|v|_inf<=30):", len(cand))
for u in cand[:5]:
    pr("   ", np.round(u, 9).tolist())
open(fn.replace("_basis.json", "_bs2_eps%g_out.txt" % eps), "w").write("\n".join(out) + "\n")
