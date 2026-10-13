import sys, time, itertools, math
import flint
from rmodels import *
n, cuts = sylvester_cuts(3)
T = int(sys.argv[1]); L = int(sys.argv[2])
A, others = slab_box(n, cuts, T, L)
p = 2**61 - 1
t0 = time.time()
print(f"H8 T={T} L={L} |A|={len(A)}", flush=True)
# incremental: dim Poly_{<=j}|_A via rank of monomials of degree <= j (exponents <= 2T)
for m in range(1, 60):
    mons = [al for al in itertools.product(range(min(2*T, m-1) + 1), repeat=7) if sum(al) < m]
    rows = [[math.prod(k[x] ** a for x, a in zip(others, al)) % p for k in A] for al in mons]
    # rank of (#mons x |A|); compute on transpose if wider
    M = flint.nmod_mat(len(rows), len(A), [v for r in rows for v in r], p)
    r = M.rank()
    print(f"  m={m}: rows={len(rows)} rank_p={r} dimV_p={len(A)-r} ({time.time()-t0:.0f}s)", flush=True)
    if r == len(A): break
