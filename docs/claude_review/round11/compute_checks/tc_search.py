"""Real tangent cone of a T=1 model, searched with the trivial strata removed (evidence + exact checks).
Since A is in {-1,0,1}^d, k_x^3 - k_x vanishes on A, so every Phi_c only involves the monomials u^alpha with
alpha <= 2 (entrywise), |alpha| = m; let mu(u) be that monomial vector.  Psi(u) = |C mu(u)|^2 / |mu(u)|^2,
C = exact integer form matrix (rows = forms, restricted to those monomials; the other coefficients are checked
to be zero).  Zeros of Psi are the TC points with mu(u) != 0.  Minimise Psi on S^{d-1} from random starts,
refine in mpmath (50 digits), and classify each zero: support (coordinates that vanish), integer relations
(LLL at 40 digits), membership in known rational subspaces.  Also checks exactly which rational subspaces
spanned by small integer vectors lie in TC."""
import sys, json, itertools, math
import numpy as np
from scipy.optimize import minimize
import mpmath as mp
import flint

fn = sys.argv[1]; starts = int(sys.argv[2]); seed = int(sys.argv[3])
F = json.load(open(fn))
mons = [tuple(a) for a in F['mons']]; rows = F['rows']; m = F['m']; d = F['d']
keep = [i for i, a in enumerate(mons) if max(a) <= 2]
for r in rows:
    assert all(r[i] == 0 for i in range(len(mons)) if i not in set(keep)), "form uses a monomial with exponent >= 3"
M2 = [mons[i] for i in keep]
C = np.array([[r[i] for i in keep] for r in rows], dtype=float)
# orthonormalise rows of C for conditioning
Qc, Rc = np.linalg.qr(C.T)
Cn = Qc.T                      # rows: orthonormal basis of the row space
E = np.array(M2, dtype=float)  # monomial exponents

def mu(u):
    return np.prod(np.power(u[None, :], E), axis=1)

def Psi(u):
    v = mu(u)
    n2 = v @ v
    if n2 == 0:
        return 1.0
    w = Cn @ v
    return (w @ w) / n2

out = []
def pr(*a):
    s = " ".join(str(x) for x in a); print(s, flush=True); out.append(s)
pr("forms %d, monomials alpha<=2 of degree %d: %d" % (len(rows), m, len(M2)))
rng = np.random.default_rng(seed)
zeros = []
vals = []
for s in range(starts):
    u0 = rng.normal(size=d); u0 /= np.linalg.norm(u0)
    r = minimize(lambda z: Psi(z / np.linalg.norm(z)), u0, method='BFGS', options=dict(gtol=1e-13, maxiter=3000))
    z = r.x / np.linalg.norm(r.x)
    vals.append(Psi(z))
    if Psi(z) < 1e-14:
        zeros.append(z)
vals.sort()
pr("Psi local minima: smallest", ["%.1e" % v for v in vals[:5]], " median %.2e" % vals[len(vals)//2], " #zeros(<1e-14):", len(zeros))

# mpmath refinement (Gauss-Newton on the forms with normalisation |u|=1 and fixed sign)
mp.mp.dps = 50
Cmp = [[mp.mpf(int(x)) for x in r] for r in [[rows[j][i] for i in keep] for j in range(len(rows))]]
def forms_mp(u):
    mv = [mp.fprod([u[t] ** a[t] for t in range(d) if a[t]]) for a in M2]
    return [mp.fsum(c * x for c, x in zip(row, mv)) for row in Cmp], mv
def refine(z):
    u = [mp.mpf(float(x)) for x in z]
    for it in range(40):
        f, mv = forms_mp(u)
        nrm = mp.sqrt(mp.fsum(x * x for x in mv))
        res = [x / nrm for x in f] + [mp.fsum(x * x for x in u) - 1]
        # numerical Jacobian
        h = mp.mpf(10) ** -25
        J = []
        for t in range(d):
            up = list(u); up[t] += h
            fp, mvp = forms_mp(up)
            nrp = mp.sqrt(mp.fsum(x * x for x in mvp))
            rp = [x / nrp for x in fp] + [mp.fsum(x * x for x in up) - 1]
            J.append([(a - b) / h for a, b in zip(rp, res)])
        Jm = mp.matrix(len(res), d)
        for i in range(len(res)):
            for t in range(d):
                Jm[i, t] = J[t][i]
        rv = mp.matrix(res)
        # least squares step via normal equations with tiny Tikhonov
        JT = Jm.T
        step = mp.lu_solve(JT * Jm + mp.eye(d) * mp.mpf(10) ** -40, JT * rv)
        u = [u[t] - step[t] for t in range(d)]
        if mp.norm(rv) < mp.mpf(10) ** -40:
            break
    f, mv = forms_mp(u)
    nrm = mp.sqrt(mp.fsum(x * x for x in mv))
    return u, mp.norm(mp.matrix(f)) / nrm

def relations_mp(u, digits=35):
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
    rel = relations_mp(u)
    rk = np.linalg.matrix_rank(np.array(rel, dtype=float)) if rel else 0
    key = (len(supp), rk)
    classes[key] = classes.get(key, 0) + 1
    if rk == 0:
        nondeg.append((float(res), [float(x) for x in u]))
pr("refined zeros classified by (support size, rank of integer-relation lattice):", classes)
pr("zeros with no integer relation (NONDEGENERATE):", len(nondeg))
for r_, u in nondeg[:5]:
    pr("   ", "%.2e" % r_, np.round(u, 8).tolist())
open(fn.replace("_forms.json", "_tcsearch_out.txt"), "w").write("\n".join(out) + "\n")
