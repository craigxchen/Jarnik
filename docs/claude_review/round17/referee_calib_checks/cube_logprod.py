#!/usr/bin/env python3
"""Referee: ln Psi_cube on C_1 by the direct product prod_x (1 + (2/3)u_x^2 t + u_x^4 t^2/18) in log domain
(all coefficients positive: no cancellation). float64 logaddexp; evidence-grade cross-check of cube_regime.py."""
import sys, numpy as np
from math import lgamma, log, cos, pi
def lnpsi(q, m, theta=0.0):
    L = np.full(m+1, -np.inf); L[0] = 0.0
    for x in range(q):
        u2 = cos(theta + 2*pi*x/q)**2
        la = log(2/3*u2) if u2 > 0 else -np.inf; lb = log(u2*u2/18) if u2 > 0 else -np.inf
        new = L.copy()
        new[1:] = np.logaddexp(new[1:], la + L[:-1])
        new[2:] = np.logaddexp(new[2:], lb + L[:-2])
        L = new
    return 2*lgamma(m+1) + L[m]
for a in sys.argv[1:]:
    q, m = map(int, a.split(':'))
    print(f'q={q} m={m}: ln Psi_cube = {lnpsi(q,m):.6f} (theta=0), {lnpsi(q,m,0.3/q):.6f} (theta=0.3/q)', flush=True)
