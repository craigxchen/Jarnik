#!/usr/bin/env python3
"""Referee: exact check of V_H = V(A,m) for small H (integer kernel basis of M^T, LLL, exact rank of result)."""
import sys, numpy as np, flint, time
from exact_full import monomials, mono_matrix
q, L, m = map(int, sys.argv[1:4]); t0=time.time()
A = np.load(f'data/A_q{q}_L{L}.npy'); M = mono_matrix(A, monomials(q, m))
MT = flint.fmpz_mat(M.T.tolist())
N, nul = MT.nullspace()            # columns of N span ker(M^T) = V(A,m)
B = flint.fmpz_mat([[int(N[i, j]) for i in range(N.nrows())] for j in range(nul)])  # rows = basis vectors
Bl = B.lll()
H = max(abs(int(Bl[i, j])) for i in range(Bl.nrows()) for j in range(Bl.ncols()))
rk = Bl.rank()
ok = (flint.fmpz_mat(M.T.tolist()) * Bl.transpose()).is_zero()
print(f'q={q} L={L} m={m}: dim V={nul}  LLL basis rank={rk}  in V: {ok}  max|c|_inf={H}  ({time.time()-t0:.0f}s)')
