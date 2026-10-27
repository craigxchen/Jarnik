"""Independent full-space (no symmetry used) check, floating point (numpy lstsq / QR), small boxes only.
Builds all of A from the orbit representatives, the full monomial evaluation matrix of degree < m (exponents <= 2),
and computes Psi_A(u) = min_P E_{mu_A}(l_u^m - P)^2 at several theta on C_1 (unnormalised u_theta), at random u
with |u|_inf = 1, plus dim V(A,m) = |A| - rank.  Also Cheb bracket [||r||^2/||r||_1, max|r|] and the exact Chebyshev
distance by LP (HiGHS) when |A| is small.
usage: python3 full_check.py q L m [lp]
"""
import sys, itertools, math, time
import numpy as np
from scipy.optimize import linprog

q, L, m = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
do_lp = len(sys.argv) > 4
reps = np.fromfile(f"data/A_q{q}_L{L}.bin", dtype=np.int8).reshape(-1, q).astype(np.int64)
pts = set()
for r in reps:
    for s in range(q): pts.add(tuple(np.roll(r, s)))
A = np.array(sorted(pts), dtype=np.float64); U = len(A)
mons = [a for deg in range(m) for a in itertools.product((0, 1, 2), repeat=q) if sum(a) == deg] if q <= 11 else None
if mons is None:
    # support-restricted enumeration
    smax = int((A != 0).sum(1).max()); mons = []
    for deg in range(m):
        for sz in range(0, min(deg, smax) + 1):
            for S in itertools.combinations(range(q), sz):
                for ex in itertools.product((1, 2), repeat=sz):
                    if sum(ex) == deg:
                        a = [0] * q
                        for x, e in zip(S, ex): a[x] = e
                        mons.append(tuple(a))
M = np.ones((U, len(mons)))
for j, a in enumerate(mons):
    for x, e in enumerate(a):
        if e: M[:, j] *= A[:, x] ** e
Qm, Rm, piv = __import__('scipy.linalg', fromlist=['qr']).qr(M, mode='economic', pivoting=True)
dR = np.abs(np.diag(Rm)); rank = int((dR > 1e-9 * dR[0]).sum()); Qm = Qm[:, :rank]
print(f"sdc{q} L={L} m={m}: |A|={U}, monomials {len(mons)}, rank {rank}, dim V = {U-rank}")


def psi(u):
    g = (A @ u) ** m
    r = g - Qm @ (Qm.T @ g)
    return (r @ r) / U, r, g


for th in np.linspace(0, 2 * np.pi / q, 5):
    u = np.cos(th + 2 * np.pi * np.arange(q) / q)
    P, r, g = psi(u)
    br = ((r @ r) / np.abs(r).sum() if np.abs(r).sum() > 0 else 0.0, np.abs(r).max())
    line = f"theta={th:.4f}: Psi={P:.12e}  E l^2m={g@g/U:.6e}  Cheb in [{br[0]:.6e}, {br[1]:.6e}]  |u|_inf={np.abs(u).max():.6f}"
    if do_lp and rank < 3000:
        # min t s.t. -t <= g - Q c <= t ; variables c (rank), t
        n = rank
        Aub = np.vstack([np.hstack([-Qm, -np.ones((U, 1))]), np.hstack([Qm, -np.ones((U, 1))])])
        bub = np.concatenate([-g, g])
        cost = np.zeros(n + 1); cost[-1] = 1
        res = linprog(cost, A_ub=Aub, b_ub=bub, bounds=[(None, None)] * n + [(0, None)], method='highs')
        line += f"  Cheb(LP)={res.fun:.8e}"
    print(line, flush=True)
rng = np.random.default_rng(1)
vals = []
for _ in range(5):
    u = rng.uniform(-1, 1, q); u /= np.abs(u).max()
    vals.append(psi(u)[0])
print("random u (|u|_inf = 1): Psi =", " ".join(f"{v:.6e}" for v in vals))
