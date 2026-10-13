# Exact integer kernel basis of V(A,m) for the H8 models, computed from scratch; LLL-reduced; saved as json.
import sys, time, json
import flint
from rmodels import *
L = int(sys.argv[1]); m = int(sys.argv[2])
n, cuts = sylvester_cuts(3)
A, others = slab_box(n, cuts, 1, L)
rows, mons = moment_matrix(A, others, m, 1)
t0 = time.time()
M = flint.fmpz_mat(rows)
r = M.rank()
print(f"L={L} m={m} |A|={len(A)} E={len(rows)} exact rank={r} dimV={len(A)-r} ({time.time()-t0:.1f}s)", flush=True)
N, nul = M.nullspace()
print("nullity", nul, f"({time.time()-t0:.1f}s)", flush=True)
# columns 0..nul-1 of N are a basis of the kernel
B = [[int(N[i, j]) for i in range(len(A))] for j in range(nul)]
# saturate: kernel lattice = Z^U cap ker. Use HNF trick: LLL on B then check primitive via the gcd of maximal minors is hard;
# simple approach: LLL reduce; then saturation via solving over Q for short vectors is not needed for Bs (span only).
Bm = flint.fmpz_mat(B)
Bl = Bm.lll()
B = [[int(Bl[i, j]) for j in range(len(A))] for i in range(Bl.nrows())]
# exact verification
Mt = M
for c in B:
    v = Mt * flint.fmpz_mat([[x] for x in c])
    assert all(int(v[i, 0]) == 0 for i in range(v.nrows()))
print("all kernel vectors verified exactly; sup norms", [max(abs(x) for x in c) for c in B])
json.dump({'A': A, 'others': others, 'basis': B, 'L': L, 'm': m}, open(f"rk_had8_L{L}_m{m}.json", "w"))
