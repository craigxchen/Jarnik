# H8 Sylvester profile, T = 2: sizes, ranges and rho (rank mod p, then exact if small) for L = 1, 2.
import sys, time
from rmodels import *
n, cuts = sylvester_cuts(3)
T = 2
for L in [int(a) for a in sys.argv[1:]]:
    t0 = time.time()
    A, others = slab_box(n, cuts, T, L)
    rngs = sorted(set(max(sum(k[x] for x in S) for k in A) - min(sum(k[x] for x in S) for k in A) for S in cuts))
    print(f"T=2 L={L}: |A|={len(A)} ranges k(S_p)={rngs} (h_A = 2*range)", flush=True)
    res = {}
    for m in range(1, 40):
        rows, _ = moment_matrix(A, others, m, T)
        r = rank_mod(rows)
        res[m] = r
        print(f"   m={m}: rank_p Poly_<m|_A = {r}, dimV_p = {len(A)-r} ({time.time()-t0:.0f}s)", flush=True)
        if r == len(A):
            break
