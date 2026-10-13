# EVIDENCE. LP-supercritical but NOT pair-separated model (as in single_tc.py): global random 0/1 profile with r columns
# on M points; sub-clusters of size n=4, all bases, slab-box T=4, w=6, m=13 (admissible iff m > LP and dim V > 0).
# Compare tangent-cone zero sets: single systems (lifted) vs the combined system, classify zeros by PSLQ relations.
import sys, itertools, random, math, numpy as np, mpmath as mp
from scipy.optimize import minimize
from engine import System
seed = int(sys.argv[1]) if len(sys.argv) > 1 else 1
M = int(sys.argv[2]) if len(sys.argv) > 2 else 6; r = int(sys.argv[3]) if len(sys.argv) > 3 else 4
random.seed(seed); rng = np.random.default_rng(seed)
n, T, w, m = 4, 4, 6, 13
cols = []
while len(cols) < r:
    a = tuple(random.randint(0, 1) for _ in range(M))
    if 0 < sum(a) < M: cols.append(a)
print(f"seed={seed} M={M} r={r} cols={cols}", flush=True)
systems = []
for J in itertools.combinations(range(M), n):
    for b in J:
        S = System(list(J), b, cols, T, w, m)
        LP = sum(S.h)/r
        if S.dimV > 0 and m > LP: systems.append(S)
print(f"admissible (J,base) systems: {len(systems)} of {math.comb(M,n)*n}", flush=True)
def tg(S, uJ):
    a = S.A @ uJ; f = a**m; nf = f @ f; g = f - S.Q @ (S.Q.T @ f); ng = g @ g
    df = (m * a**(m-1))[:, None] * S.A
    return ng/nf, (2*(g @ df)*nf - ng*2*(f @ df))/nf**2
def obj(x, syss):
    u = np.append(x, 0.0); u -= u.mean(); nu = np.linalg.norm(u); u /= nu
    val = 0.0; gu = np.zeros(M)
    for S in syss:
        t, g = tg(S, u[S.J]); val += t; gu[S.J] += g
    P = np.eye(M) - 1.0/M; G = (np.eye(M) - np.outer(u, u)) @ gu / nu
    return val/len(syss), (P @ G)[:-1]/len(syss)
def rel(u):
    # integer relation among the differences u_x - u_last  <=>  v in Z^n_0 with <v,u> = 0 (v_last = -sum)
    mp.mp.dps = 30
    d = [mp.mpf(float(x - u[-1])) for x in u[:-1]]
    v = mp.pslq(d, tol=1e-8, maxcoeff=30, maxsteps=10**5)
    return None if v is None else list(v) + [-sum(v)]
def zeros(syss, starts, label):
    out = []
    for s in range(starts):
        r_ = minimize(obj, rng.normal(size=M-1), args=(syss,), jac=True, method='BFGS', options={'gtol': 1e-16, 'maxiter': 5000})
        u = np.append(r_.x, 0.0); u -= u.mean(); u /= np.linalg.norm(u)
        if u[np.argmax(np.abs(u))] < 0: u = -u
        out.append((r_.fun, u))
    out.sort(key=lambda t: t[0])
    zs = []
    for v, u in out:
        if v < 1e-20 and not any(np.linalg.norm(u - z) < 1e-6 for z in zs): zs.append(u)
    rels = [rel(z) for z in zs]
    print(f"[{label}] min residual {out[0][0]:.1e}, median {out[len(out)//2][0]:.1e}; distinct zeros {len(zs)}; "
          f"with an integer relation: {sum(1 for x in rels if x is not None)}", flush=True)
    for z, rr in list(zip(zs, rels))[:8]:
        print(f"    u/max={np.round(z/np.abs(z).max(), 4)} relation={rr}", flush=True)
    return zs, rels
# single systems: zeros in R^M where only u_J matters -> restrict to R^J by optimising over all of R^M is degenerate;
# instead lift: minimise single tau over u_J and report the J-projection classification
nondeg_single = 0; tot_single = 0
for S in systems[:int(sys.argv[4]) if len(sys.argv) > 4 else 12]:
    best = []
    for s in range(20):
        def o(x):
            u = np.append(x, 0.0); u -= u.mean(); nu = np.linalg.norm(u); u /= nu
            t, g = tg(S, u); P = np.eye(n) - 1.0/n; G = (np.eye(n) - np.outer(u, u)) @ g / nu
            return t, (P @ G)[:-1]
        r_ = minimize(o, rng.normal(size=n-1), jac=True, method='BFGS', options={'gtol': 1e-16, 'maxiter': 5000})
        u = np.append(r_.x, 0.0); u -= u.mean(); u /= np.linalg.norm(u)
        if u[np.argmax(np.abs(u))] < 0: u = -u
        best.append((r_.fun, u))
    zs = []
    for v, u in sorted(best, key=lambda t: t[0]):
        if v < 1e-20 and not any(np.linalg.norm(u - z) < 1e-6 for z in zs): zs.append(u)
    rr = [rel(z) for z in zs]
    tot_single += len(zs); nondeg_single += sum(1 for x in rr if x is None)
    print(f"  single J={S.J} base={S.base} dimV={S.dimV}: zeros {len(zs)}, without relation {sum(1 for x in rr if x is None)}", flush=True)
print(f"single systems: total TC zeros {tot_single}, without integer relation (|coef|<=30) {nondeg_single}", flush=True)
zeros(systems, int(sys.argv[5]) if len(sys.argv) > 5 else 40, "combined")

# ---- margin-constrained minimisation (added): inf of the TC residual over eps-nondegenerate directions ----
# Region R_eps = {u in S(R^K_0) : |<v,u>| >= eps |v|_1 for all 0 != v in Z^K_0, |v|_inf <= 2}.  Minimise
# tau(u) + 1e3 * sum_v max(0, eps - |<v,u>|/|v|_1)^2 from many random starts; report min tau inside R_eps.
def vset(K, vmax=2):
    out = []
    for v in itertools.product(range(-vmax, vmax+1), repeat=K-1):
        if any(v): vv = np.array(list(v) + [-sum(v)], dtype=float); out.append(vv/np.abs(vv).sum())
    return np.array(out)
def margin_min(syss, K, eps, starts, local_index):
    V = vset(K)
    def ob(x):
        u = np.append(x, 0.0); u -= u.mean(); nu = np.linalg.norm(u); u /= nu
        val = 0.0; gu = np.zeros(K)
        for S in syss:
            idx = [local_index[j] for j in S.J]
            t, g = tg(S, u[idx]); val += t; gu[idx] += g
        val /= len(syss); gu /= len(syss)
        pr = V @ u; viol = np.maximum(0, eps - np.abs(pr))
        val += 1e3*np.sum(viol**2); gu += 1e3*(-2*viol*np.sign(pr)) @ V
        P = np.eye(K) - 1.0/K; G = (np.eye(K) - np.outer(u, u)) @ gu / nu
        return val, (P @ G)[:-1]
    best = np.inf
    for s in range(starts):
        r_ = minimize(ob, rng.normal(size=K-1), jac=True, method='BFGS', options={'gtol': 1e-14, 'maxiter': 4000})
        u = np.append(r_.x, 0.0); u -= u.mean(); u /= np.linalg.norm(u)
        if np.min(np.abs(V @ u)) >= 0.999*eps:
            t = np.mean([tg(S, u[[local_index[j] for j in S.J]])[0] for S in syss]); best = min(best, t)
    return best
for eps in [0.05, 0.02]:
    vals = []
    for S in systems[:8]:
        li = {j: i for i, j in enumerate(S.J)}
        vals.append(margin_min([S], n, eps, 25, li))
    print(f"eps={eps}: single systems (8 sampled): min tau over R_eps = {[f'{v:.1e}' for v in vals]}", flush=True)
    li = {j: j for j in range(M)}
    print(f"eps={eps}: combined ({len(systems)} systems): min tau over R_eps = {margin_min(systems, M, eps, 25, li):.2e}", flush=True)
