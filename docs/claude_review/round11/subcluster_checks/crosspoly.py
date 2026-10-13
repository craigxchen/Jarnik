# EXACT (mod p ranks). Cross-polytope profiles: points {0, +-e_i} in Z^rho (M = 2rho+1, rank rho),
# columns <eps, lambda> for eps in {+-1}^rho modulo sign (r = 2^(rho-1)), e_p = 2, omega_p = 2^(1-rho).
# Pair separation is tight: ||lambda_x - lambda_y|| = E_eps |<eps, v>| >= 1.  We compute, for sub-clusters J of size n
# (all bases), the pushforward dims, the combined dim, and dim U = V(S_slab, m) of the universal slab system.
import sys, itertools
from pushforward2 import Profile, run
rho = int(sys.argv[1])
pts = [tuple([0]*rho)] + [tuple((s if j == i else 0) for j in range(rho)) for i in range(rho) for s in (1, -1)]
dirs = [eps for eps in itertools.product([1, -1], repeat=rho) if eps[0] == 1]
P = Profile(pts, dirs, [1/len(dirs)]*len(dirs))
pn = [P.norm(tuple(a-b for a, b in zip(x, y))) for x, y in itertools.combinations(pts, 2)]
print(f"rho={rho}: M={len(pts)} r={len(dirs)} e={P.e}; min pair norm {min(pn):.3f}; "
      f"min norm over Lambda=Z^rho (box 3): {min(P.norm(l) for l in itertools.product(range(-3,4), repeat=rho) if any(l)):.3f}")
for spec in sys.argv[2:]:
    n, T, w, m = map(int, spec.split(','))
    run(P, n, T, w, m, f"cross-polytope rho={rho}", want_float=False)
