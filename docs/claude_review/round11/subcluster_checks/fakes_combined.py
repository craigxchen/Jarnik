# EVIDENCE (floating point). Full-rank (fake-type) global profiles on M points; sub-clusters J of size n with ALL base
# points; slab-box per (J,base) from the restricted columns (median cut as in the round-10 evidence scripts; these
# small models are NOT LP-supercritical: they probe the geometry of the base loci only),
# m = largest with 4*binom(m+n-2,n-1) <= |A|.
# (1) fake positions delta ~ U[0,Delta]^M and random positions: certified off every single Z_{P_J}
#     (|leading Taylor term| > Taylor remainder bound), counted per position;
# (2) tangent cones: combined TC residual  mean_J tau_J(u_J)  minimised on the unit sphere of R^M_0 (all (J,base));
#     compared with the minima of single-J tau over u in R^M_0 (u_J only matters). Zeros classified by PSLQ.
import sys, itertools, math, numpy as np, mpmath as mp
from scipy.optimize import minimize
src = open('/home/user/Jarnik/docs/claude_review/round10/referee_determinant_checks/ref_fakes.py').read().split('\nout=[]')[0]
exec(src)            # defines paley_I, skew_double_core, round5_profile (read-only reuse of referee code)
from engine import poly_basis
rng = np.random.default_rng(11)

def build(Aexp, n, T):
    M, r = Aexp.shape
    beta_all = (2*Aexp - 1).T.astype(float)       # r x M, e_p = 1
    systems = []
    for J in itertools.combinations(range(M), n):
        for b in J:
            bi = J.index(b); others = [i for i in range(n) if i != bi]
            box = np.array(list(itertools.product(range(-T, T+1), repeat=n-1)))
            K = np.zeros((len(box), n), dtype=int); K[:, others] = box; K[:, bi] = -box.sum(1)
            X = K @ beta_all[:, list(J)].T
            mx = np.abs(X).max(axis=1); w = np.sort(mx)[len(mx)//2]
            A = K[mx <= w]
            m = 1
            while 4*math.comb(m+1+n-2, n-1) <= len(A): m += 1
            Q, rk, gap = poly_basis([tuple(int(t) for t in k) for k in A], bi, m, T)
            systems.append(dict(J=J, b=b, A=A.astype(float), m=m, Q=Q, dimV=len(A)-rk))
    return systems

def tau_and_grad(S, uJ):
    a = S['A'] @ uJ; m = S['m']
    f = a**m; nf = f @ f
    g = f - S['Q'] @ (S['Q'].T @ f); ng = g @ g
    df = (m * a**(m-1))[:, None] * S['A']            # d f / d uJ
    dng = 2 * (g @ df); dnf = 2 * (f @ df)
    return ng/nf, (dng*nf - ng*dnf)/nf**2

def combined_obj(x, systems, M):
    u = np.append(x, 0.0); u = u - u.mean(); nu = np.linalg.norm(u); u = u/nu
    val = 0.0; grad_u = np.zeros(M)
    for S in systems:
        J = list(S['J']); t, gr = tau_and_grad(S, u[J]); val += t; grad_u[J] += gr
    val /= len(systems); grad_u /= len(systems)
    # chain rule through centring and normalisation
    P = np.eye(M) - 1.0/M; G = (np.eye(M) - np.outer(u, u)) @ grad_u / nu
    gx = (P @ G)[:-1]
    return val, gx

def certify_off(S, dJ):
    # sum_k c_k e^{i<k,d>} != 0 for some c in V(A,m): use V = complement of Q; c = projection of the degree-m Taylor vector
    a = S['A'] @ dJ; m = S['m']
    f = a**m; c = f - S['Q'] @ (S['Q'].T @ f)              # c in V(A,m) (float); Taylor: sum c_k a_k^j = 0 for j<m
    lead = abs(c @ (a**m))/math.factorial(m)
    rem = np.abs(c) @ (np.abs(a)**(m+1))/math.factorial(m+1)
    return lead > 1.01*rem, lead, rem

def relation(u):
    # integer relation among differences u_x - u_last (i.e. v in Z^M_0 with <v,u> = 0)
    mp.mp.dps = 30
    v = mp.pslq([mp.mpf(float(x - u[-1])) for x in u[:-1]], tol=1e-9, maxcoeff=50, maxsteps=10**5)
    return None if v is None else list(v) + [-sum(v)]

profiles = [("skew double core q=7", skew_double_core(7)), ("round-5 Paley-I b=8 q=7", round5_profile(7))]
which = int(sys.argv[1]) if len(sys.argv) > 1 else 0
n = int(sys.argv[2]) if len(sys.argv) > 2 else 4; T = int(sys.argv[3]) if len(sys.argv) > 3 else 3
name, Aexp = profiles[which]
M, r = Aexp.shape
print(f"{name}: M={M} r={r}; sub-clusters n={n}, T={T}, all bases", flush=True)
systems = build(Aexp, n, T)
print(f"  systems: {len(systems)}; m values {sorted(set(S['m'] for S in systems))}; |A| {min(len(S['A']) for S in systems)}..{max(len(S['A']) for S in systems)}; "
      f"dim V {min(S['dimV'] for S in systems)}..{max(S['dimV'] for S in systems)}", flush=True)
# (1) certification at fake-like and random positions
for lab, gen in [("fake-type U[0,1e-4]^M", lambda: rng.uniform(0, 1e-4, M)), ("random N(0,1e-3)", lambda: rng.normal(0, 1e-3, M))]:
    allsys = 0; per = []
    for t in range(50):
        d = gen(); d -= d.mean()
        ok = [certify_off(S, d[list(S['J'])])[0] for S in systems]
        per.append(sum(ok)); allsys += all(ok)
    print(f"  (1) {lab}: positions certified off EVERY single Z_PJ: {allsys}/50; min #systems certifying per position {min(per)}/{len(systems)}", flush=True)
# (2) combined tangent cone
res = []
for s in range(int(sys.argv[4]) if len(sys.argv) > 4 else 40):
    x0 = rng.normal(size=M-1)
    r_ = minimize(combined_obj, x0, args=(systems, M), jac=True, method='BFGS', options={'gtol': 1e-14, 'maxiter': 3000})
    u = np.append(r_.x, 0.0); u -= u.mean(); u /= np.linalg.norm(u)
    res.append((r_.fun, u))
res.sort(key=lambda t: t[0])
print(f"  (2) combined TC: smallest residuals {[f'{v:.1e}' for v, _ in res[:8]]}", flush=True)
for v, u in res[:6]:
    rel = relation(u) if v < 1e-12 else None
    print(f"     resid {v:.2e}  u/max = {np.round(u/np.abs(u).max(), 4)}  relation {rel}", flush=True)
