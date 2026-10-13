"""Search the real base locus of a T = 1 model near e and EXPLAIN every zero by an exactly verified subtorus.
System: bs_search2 formulation (tau = tan(u/2), relative residuals, Levenberg-Marquardt on spheres |tau| = eps).
For each numerical zero u: integer relations by LLL (|v|_inf <= 30, tolerance 1e-7 eps), H = their common
kernel (rational subspace).  Exact test: exp(i H) lies in Bs(A,m)  iff  for every basis c of V(A,m) the
push-forward of c under k -> (<k,g_1>,...,<k,g_j>) (g = integer basis of H) is zero (distinct characters on the
image are linearly independent).  A zero is EXPLAINED if it lies (to 1e-8 relative) on a verified subtorus.
Evidence only: the search is by random starts."""
import sys, json, itertools, math
import numpy as np
from scipy.optimize import least_squares
import flint

def build(fn):
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
    Nex = np.array([[int(np.dot(cj, col)) for col in cols] for cj in Bo], dtype=float)
    E = np.array(alphas); n1 = (E == 1).sum(1); n2 = (E == 2).sum(1)
    coef = (2.0 ** (n1 + n2)) * ((-1.0) ** n2) * np.where(n1 % 4 < 2, 1.0, -1.0)
    par = n1 % 2; supp = E > 0
    def resid(t):
        t2 = 1 + t * t
        mon = np.prod(np.where(supp, t[None, :] ** E, 1.0), axis=1)
        rest = np.prod(np.where(supp, 1.0, t2[None, :]), axis=1)
        T_ = coef * mon * rest
        out = []
        for p in (0, 1):
            mask = par == p
            num = Nex[:, mask] @ T_[mask]
            den = np.abs(Nex[:, mask]) @ np.abs(T_[mask]) + 1e-300
            out.append(num / den)
        return np.concatenate(out)
    return D, A, basis, m, d, resid

def relations(u, tol, bound=30):
    d = len(u)
    K = 10 ** 11
    rows = [[1 if j == i else 0 for j in range(d)] + [int(round(K * u[i]))] for i in range(d)]
    R = flint.fmpz_mat(rows).lll()
    rel = []
    for i in range(R.nrows()):
        v = [int(R[i, j]) for j in range(d)]
        if any(v) and max(abs(x) for x in v) <= bound and abs(np.dot(v, u)) < tol:
            rel.append(v)
    return rel

def subspace_basis(rel, d):
    if not rel:
        return [[1 if j == i else 0 for j in range(d)] for i in range(d)]
    M = flint.fmpz_mat(rel)
    X, nul = M.nullspace()
    G = flint.fmpz_mat([[int(X[i, j]) for i in range(d)] for j in range(nul)])
    if nul == 0:
        return []
    G = G.lll()
    return [[int(G[i, j]) for j in range(d)] for i in range(G.nrows()) if any(G[i, j] for j in range(d))]

def subtorus_in_Bs(A, basis, g):
    if not g:
        return True
    ys = [tuple(sum(k[t] * gg[t] for t in range(len(k))) for gg in g) for k in A]
    for c in basis:
        acc = {}
        for y, ck in zip(ys, c):
            if ck:
                acc[y] = acc.get(y, 0) + ck
        if any(v != 0 for v in acc.values()):
            return False
    return True

if __name__ == "__main__":
    fn = sys.argv[1]; starts = int(sys.argv[2]); eps_list = [float(x) for x in sys.argv[3].split(',')]
    seed = int(sys.argv[4])
    D, A, basis, m, d, resid = build(fn)
    rng = np.random.default_rng(seed)
    out = []
    def pr(*a):
        s = " ".join(str(x) for x in a); print(s, flush=True); out.append(s)
    verified = {}
    for eps in eps_list:
        stats = dict(nonzero=0, explained=0, unexplained=0)
        unexpl = []
        for s in range(starts):
            z0 = rng.normal(size=d); z0 *= eps / np.linalg.norm(z0)
            f = lambda z: np.concatenate([resid(z), [(z @ z - eps * eps) / eps ** 2]])
            r = least_squares(f, z0, method='lm', xtol=1e-15, ftol=1e-15, gtol=1e-15, max_nfev=20000)
            z = r.x
            if np.linalg.norm(resid(z)) > 1e-10:
                stats['nonzero'] += 1
                continue
            u = 2 * np.arctan(z)
            rel = relations(u, 1e-7 * eps)
            g = subspace_basis(rel, d)
            key = tuple(map(tuple, sorted(g)))
            if key not in verified:
                verified[key] = subtorus_in_Bs(A, basis, g) if len(g) < d else False
            ok = verified[key]
            # distance of u to span(g)
            if g:
                Gm = np.array(g, dtype=float).T
                proj = Gm @ np.linalg.lstsq(Gm, u, rcond=None)[0]
                dist = np.linalg.norm(u - proj) / eps
            else:
                dist = 1.0
            if ok and dist < 1e-8:
                stats['explained'] += 1
            else:
                stats['unexplained'] += 1
                unexpl.append((dist, len(g), np.round(u / eps, 6).tolist()))
        pr("eps=%g: %s" % (eps, stats))
        for x in unexpl[:6]:
            pr("   unexplained zero: dist=%.1e dimH=%d u/eps=%s" % x)
    pr("verified subtori found (dimension: count):", {})
    dims = {}
    for key, ok in verified.items():
        if ok:
            dims[len(key)] = dims.get(len(key), 0) + 1
    pr("   ", dims)
    pr("candidate subspaces REJECTED by the exact push-forward test:", sum(1 for ok in verified.values() if not ok))
    for key, ok in verified.items():
        if ok:
            out.append("   verified subtorus basis: %s" % (list(key),))
    open(fn.replace("_basis.json", "_explain_out.txt"), "w").write("\n".join(out) + "\n")
