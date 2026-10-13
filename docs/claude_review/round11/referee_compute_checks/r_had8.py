import sys, time
from rmodels import *
n, cuts = sylvester_cuts(3)
out = []
for L in [1, 2, 3]:
    t0 = time.time()
    A, others = slab_box(n, cuts, 1, L)
    res = rho_and_dims(A, others, 1, 16)
    ranges = []
    for S in cuts:
        vals = [sum(k[x] for x in S) for k in A]
        ranges.append(max(vals) - min(vals))
    rho = max(m for m, r in res.items() if r < len(A)) if any(r < len(A) for r in res.values()) else 0
    line = f"H8 T=1 L={L}: |A|={len(A)} ranges k(S_p)={sorted(set(ranges))} max h_A=2*range={2*max(ranges)} rank Poly_<m (mod p)={res} rho={rho}"
    line += "  dimV(A,m)=" + str({m: len(A) - r for m, r in res.items()}) + f" ({time.time()-t0:.1f}s)"
    print(line, flush=True)
