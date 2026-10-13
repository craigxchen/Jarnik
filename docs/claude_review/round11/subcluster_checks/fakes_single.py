# EVIDENCE. Single-(J,base) tangent cones for the same systems as fakes_combined.py: minimise tau_J over unit u in R^J_0.
import sys, numpy as np
exec(open('fakes_combined.py').read().split("profiles = [")[0])
profiles = [("skew double core q=7", skew_double_core(7)), ("round-5 Paley-I b=8 q=7", round5_profile(7))]
which = int(sys.argv[1]); n = int(sys.argv[2]); T = int(sys.argv[3]); nJ = int(sys.argv[4])
name, Aexp = profiles[which]; M = Aexp.shape[0]
systems = build(Aexp, n, T)
pick = rng.choice(len(systems), nJ, replace=False)
zero_count = 0; deg = 0; nondeg = 0
for i in pick:
    S = systems[i]
    def obj(x):
        u = np.append(x, 0.0); u -= u.mean(); nu = np.linalg.norm(u); u /= nu
        t, g = tau_and_grad(S, u)
        P = np.eye(n) - 1.0/n; G = (np.eye(n) - np.outer(u, u)) @ g / nu
        return t, (P @ G)[:-1]
    best = []
    for s in range(30):
        r_ = minimize(obj, rng.normal(size=n-1), jac=True, method='BFGS', options={'gtol': 1e-15, 'maxiter': 3000})
        u = np.append(r_.x, 0.0); u -= u.mean(); u /= np.linalg.norm(u)
        best.append((r_.fun, u))
    best.sort(key=lambda t: t[0])
    zs = []
    for v, u in best:
        if v < 1e-20:
            if u[np.argmax(np.abs(u))] < 0: u = -u
            if not any(np.linalg.norm(u - z) < 1e-6 for z in zs): zs.append(u)
    rels = [relation(z) for z in zs]
    zero_count += len(zs); deg += sum(1 for r in rels if r is not None); nondeg += sum(1 for r in rels if r is None)
    print(f"  J={S['J']} base={S['b']} m={S['m']}: min tau {best[0][0]:.2e}; zeros {len(zs)}: " +
          "; ".join(f"{np.round(z/np.abs(z).max(),3)} rel {r}" for z, r in zip(zs, rels)), flush=True)
print(f"{name}: sampled {nJ} single systems; TC zeros found {zero_count} (degenerate {deg}, no relation found {nondeg})")
