"""Exact basis of V(A,13) for the undecided n=5 instance (nullspace over Q + LLL), saved for the cone pipeline."""
import json
import flint
from core import build_A, moment_matrix
prof = ((0, 2), (0, 3), (1,), (2,), (2, 4), (3, 4), (4,))
A = build_A(5, 3, [(S, 3) for S in prof])
M = moment_matrix(A, 13)
X, nul = M.nullspace()
U = len(A)
B = flint.fmpz_mat([[X[i, j] for i in range(U)] for j in range(nul)]).lll()
basis = [[int(B[i, j]) for j in range(U)] for i in range(B.nrows())]
for c in basis:
    v = M * flint.fmpz_mat([[x] for x in c]); assert all(v[i, 0] == 0 for i in range(v.nrows()))
print("|A| =", U, "dim V(A,13) =", len(basis), "sup =", max(abs(x) for c in basis for x in c))
json.dump(dict(A=A, basis=basis, T=3, L=3, m=13, prof=prof), open("n5hard_T3_L3_m13_basis.json", "w"))
