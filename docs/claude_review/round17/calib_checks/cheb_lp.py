"""Chebyshev distance Cheb_A(u_theta) = d^m kappa_V(u_theta) on the circle C_1 by LP (HiGHS via scipy.linprog), full A,
on a theta grid of the fundamental domain [0, pi/q] (theta -> -theta is x -> -x, theta -> theta + 2pi/q is a translation,
both box automorphisms).  Floating point: EVIDENCE.  Reports kappa_V(hat u) = Cheb/(|u|_inf^m d^m), the ratio inf/sup
(barrier Thm 3.2 predicts >= (m+1)^{-1/2}), the L2 lower bound sqrt(Psi)/|u|_inf^m, and the dual measure's l1 norm.
usage: python3 cheb_lp.py q L m [ngrid]
"""
import sys, itertools, math, time
import numpy as np
import scipy.linalg as sla
from scipy.optimize import linprog
from scipy.sparse import csr_matrix, hstack as sphstack, vstack as spvstack

q, L, m = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
ng = int(sys.argv[4]) if len(sys.argv) > 4 else 9
t0 = time.time()
reps = np.fromfile(f"data/A_q{q}_L{L}.bin", dtype=np.int8).reshape(-1, q).astype(np.int64)
pts = set()
for r in reps:
    for s in range(q): pts.add(tuple(np.roll(r, s)))
A = np.array(sorted(pts), dtype=np.float64); U = len(A)
smax = int((A != 0).sum(1).max()); mons = []
for deg in range(m):
    for sz in range(0, min(deg, smax) + 1):
        for S in itertools.combinations(range(q), sz):
            for ex in itertools.product((1, 2), repeat=sz):
                if sum(ex) == deg:
                    a = [0] * q
                    for x, e in zip(S, ex): a[x] = e
                    mons.append(a)
M = np.ones((U, len(mons)))
for j, a in enumerate(mons):
    for x, e in enumerate(a):
        if e: M[:, j] *= A[:, x] ** e
_, Rm, piv = sla.qr(M, mode='economic', pivoting=True)
dR = np.abs(np.diag(Rm)); rank = int((dR > 1e-9 * dR[0]).sum())
Ms = M[:, piv[:rank]]
Qo, _ = np.linalg.qr(Ms)
print(f"sdc{q} L={L} m={m}: |A|={U}, rank Poly_<m = {rank}, dim V = {U-rank}  ({time.time()-t0:.0f}s)", flush=True)
Ssp = csr_matrix(Ms)
rows = []
for th in np.linspace(0, np.pi / q, ng):
    u = np.cos(th + 2 * np.pi * np.arange(q) / q); ui = np.abs(u).max()
    g = (A @ u) ** m
    r = g - Qo @ (Qo.T @ g); psi = r @ r / U
    one = csr_matrix(np.ones((U, 1)))
    Aub = spvstack([sphstack([-Ssp, -one]), sphstack([Ssp, -one])]).tocsr()
    bub = np.concatenate([-g, g])
    cost = np.zeros(rank + 1); cost[-1] = 1
    res = linprog(cost, A_ub=Aub, b_ub=bub, bounds=[(None, None)] * rank + [(0, None)], method='highs')
    ch = res.fun
    lam = res.ineqlin.marginals if hasattr(res, 'ineqlin') else None
    dual_l1 = float(np.abs(lam).sum()) if lam is not None else float('nan')
    kap = ch / (ui ** m * q ** m)
    rows.append((th, ch, ui, kap, math.sqrt(psi) / ui ** m))
    print(f"theta={th:.5f}: Cheb={ch:.8e} |u|_inf={ui:.6f} kappa_V(hat u)={kap:.6e}  sqrt(Psi(hat u))={math.sqrt(psi)/ui**m:.6e}"
          f"  Cheb/sqrt(Psi)={ch/math.sqrt(psi):.4f}  dual l1={dual_l1:.4f}", flush=True)
ks = [r[3] for r in rows]
print(f"inf/sup kappa_V(hat u) over grid = {min(ks)/max(ks):.5f}  (Thm 3.2 bound (m+1)^-1/2 = {(m+1)**-0.5:.5f});"
      f" inf kappa_V = {min(ks):.6e} at theta={rows[int(np.argmin(ks))][0]:.5f}; ln inf kappa_V = {math.log(min(ks)):.4f}  ({time.time()-t0:.0f}s)")
