"""Variant support for the Hadamard profile at n = 8 (README 'freedom'): the symmetric box |k_x| <= 1 for ALL
eight points (also the base point, k_8 = -sum k'), slabs |k(S_p)| <= L.  Exact rho, LP bound on zeta (with the
extra functional), exact basis of V(A,m), saved for the cone pipeline."""
import sys, json
import flint
from core import build_A, moment_matrix, width, siegel_H
from bigger import hadamard8
from lpzeta import lp_bound
from zono import functionals
L = int(sys.argv[1])
prof = hadamard8()
A = [k for k in build_A(8, 1, [(S, L) for S in prof]) if abs(sum(k)) <= 1]
U = len(A); j = 0
while moment_matrix(A, j + 1).rank() < U:
    j += 1
rho = j
lams = functionals(8, 1, prof, L) + [(tuple([1] * 7), 1)]
lb, _ = lp_bound(None, 1, None, None, d=7, lams=lams)
print("had8 symmetric box T=1 L=%d: |A|=%d rho=%d LP zeta<=%s ranges=%s" % (L, U, rho, lb, [width(A, S) for S in prof]))
m = rho
if not lb < m:
    print("not faithful"); sys.exit(0)
M = moment_matrix(A, m)
X, nul = M.nullspace()
B = flint.fmpz_mat([[X[i, jj] for i in range(U)] for jj in range(nul)]).lll()
basis = [[int(B[i, jj]) for jj in range(U)] for i in range(B.nrows())]
for c in basis:
    v = M * flint.fmpz_mat([[x] for x in c]); assert all(v[i, 0] == 0 for i in range(v.nrows()))
print("m=%d dimV=%d sup=%d H=%.3g" % (m, len(basis), max(abs(x) for c in basis for x in c), siegel_H(A, 1, m)))
json.dump(dict(A=A, basis=basis, T=1, L=L, m=m, prof=prof), open("had8sym_T1_L%d_m%d_basis.json" % (L, m), "w"))
