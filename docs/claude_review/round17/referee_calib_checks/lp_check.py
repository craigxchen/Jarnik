#!/usr/bin/env python3
"""Referee: kappa_V(hat u_theta) = Cheb_A(hat u_theta)/d^m by LP (HiGHS), independent code. Evidence (floating LP)."""
import sys, numpy as np, itertools
from scipy.optimize import linprog
from math import cos, pi
sys.path.insert(0, '.')
from exact_full import monomials, mono_matrix
q, L, m = map(int, sys.argv[1:4]); nth = int(sys.argv[4]) if len(sys.argv) > 4 else 5
A = np.load(f'data/A_q{q}_L{L}.npy')
mons = monomials(q, m); M = mono_matrix(A, mons).astype(float)
# column basis (float QR with pivoting)
import scipy.linalg as sl
Q_, R_, P_ = sl.qr(M, mode='economic', pivoting=True)
rk = int((np.abs(np.diag(R_)) > 1e-8*abs(R_[0,0])).sum()); B = M[:, P_[:rk]]
U = A.shape[0]; vals = []
for th in np.linspace(0, pi/(2*q), nth):
    u = np.array([cos(th + 2*pi*x/q) for x in range(q)]); u /= np.abs(u).max()
    g = (A @ u)**m
    # min t s.t. -t <= g - B c <= t ; vars (c, t)
    c_obj = np.zeros(rk+1); c_obj[-1] = 1
    Aub = np.vstack([np.hstack([-B, -np.ones((U,1))]), np.hstack([B, -np.ones((U,1))])])
    bub = np.concatenate([-g, g])
    r = linprog(c_obj, A_ub=Aub, b_ub=bub, bounds=[(None,None)]*rk+[(0,None)], method='highs')
    vals.append(r.fun / q**m)
    print(f'theta={th:.5f}: Cheb={r.fun:.8g}  kappa_V={r.fun/q**m:.6e}', flush=True)
print(f'q={q} L={L} m={m}: inf kappa_V={min(vals):.4e}  inf/sup={min(vals)/max(vals):.4f}')
