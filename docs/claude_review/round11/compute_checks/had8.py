"""Faithful model: Sylvester-Hadamard profile at n = 8 (d = 7), box T = 1, slab half-width L, order m.
(1) exact rho (fmpz rank of moment matrices), exact integer kernel basis of V(A,m) (nullspace + LLL),
    sup-norm vs Siegel H.  Saves basis to had8_T{T}_L{L}_m{m}_basis.txt.
"""
import sys, time, json
import flint
from core import build_A, moment_matrix, width, siegel_H
from bigger import hadamard8

T = int(sys.argv[1]); L = int(sys.argv[2]); m = int(sys.argv[3])
n = 8; prof = hadamard8()
A = build_A(n, T, [(S, L) for S in prof])
out = []
def pr(*a):
    s = " ".join(str(x) for x in a); print(s, flush=True); out.append(s)
pr("had8 T=%d L=%d m=%d |A|=%d profile=%s" % (T, L, m, len(A), prof))
pr("ranges of k(S_p):", [width(A, S) for S in prof])
t0 = time.time()
M = moment_matrix(A, m)
r = M.rank()
M1 = moment_matrix(A, m + 1)
r1 = M1.rank()
pr("exact rank Poly_{<m}|A =", r, " Poly_{<m+1}|A =", r1, " |A| =", len(A), " (%.1fs)" % (time.time() - t0))
pr("dim V(A,m) =", len(A) - r, "; dim V(A,m+1) =", len(A) - r1, "; rho >= m:", r < len(A), "; rho = m:", r1 == len(A))
t0 = time.time()
X, nul = M.nullspace()
U = len(A)
B = flint.fmpz_mat([[X[i, j] for i in range(U)] for j in range(nul)])
pr("nullspace computed (%.1fs), raw sup %d" % (time.time() - t0, max(abs(int(x)) for x in B.entries())))
t0 = time.time()
B = B.lll()
sup = max(abs(int(x)) for x in B.entries())
pr("LLL (%.1fs): basis of a full-rank sublattice of V(A,m) cap Z^A, sup-norm %d, H = %.3g, V_H = V(A,m): %s" % (
    time.time() - t0, sup, siegel_H(A, T, m), sup <= siegel_H(A, T, m)))
basis = [[int(B[i, j]) for j in range(U)] for i in range(B.nrows())]
# verify exactly
for c in basis:
    v = M * flint.fmpz_mat([[x] for x in c])
    assert all(v[i, 0] == 0 for i in range(v.nrows()))
pr("all", len(basis), "basis vectors verified exactly in V(A,m)")
json.dump(dict(A=A, basis=basis, T=T, L=L, m=m, prof=prof), open("had8_T%d_L%d_m%d_basis.json" % (T, L, m), "w"))
open("had8_T%d_L%d_m%d_out.txt" % (T, L, m), "w").write("\n".join(out) + "\n")
