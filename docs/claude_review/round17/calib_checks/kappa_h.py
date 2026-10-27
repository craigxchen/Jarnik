"""kappa_H on the circle C_1 for small toy boxes (feasible only for small |A|, dim V).
Exact part: integer basis of V(A,m) = ker (moment matrix of degree < m)^T over Z (flint fmpz_mat nullspace, primitive
rows), LLL-reduced (flint).  H_red = max |c|_inf over the reduced basis: for H >= H_red, V_H = V(A,m) (exact), so
kappa_H = kappa_V and Psi_{A,H} = Psi_A.  For smaller H the span of the reduced vectors with |c|_inf <= H is a subspace
of V_H, which gives the lower bound  d^m kappa_H(u) >= Psi_{A,span}(u)^{1/2}  (reduction §4, same proof as Thm 2.1(2);
the projection is floating point: EVIDENCE for the numbers, exact for the spans).
usage: python3 kappa_h.py q L m
"""
import sys, itertools, math, time
import numpy as np
import flint

q, L, m = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
t0 = time.time()
reps = np.fromfile(f"data/A_q{q}_L{L}.bin", dtype=np.int8).reshape(-1, q).astype(np.int64)
pts = set()
for r in reps:
    for s in range(q): pts.add(tuple(np.roll(r, s)))
A = np.array(sorted(pts), dtype=np.int64); U = len(A)
smax = int((A != 0).sum(1).max()); mons = []
for deg in range(m):
    for sz in range(0, min(deg, smax) + 1):
        for S in itertools.combinations(range(q), sz):
            for ex in itertools.product((1, 2), repeat=sz):
                if sum(ex) == deg:
                    a = [0] * q
                    for x, e in zip(S, ex): a[x] = e
                    mons.append(a)
M = np.ones((U, len(mons)), dtype=np.int64)
for j, a in enumerate(mons):
    for x, e in enumerate(a):
        if e: M[:, j] *= A[:, x] ** e
# independent monomial columns (rank mod a large prime = rank over Q generically; we use exact fmpz rank below)
Mt = flint.fmpz_mat(M.T.tolist())
rk = Mt.rank()
print(f"sdc{q} L={L} m={m}: |A|={U}, monomials {len(mons)}, rank {rk}, dim V = {U-rk}  ({time.time()-t0:.0f}s)", flush=True)
X, nul = Mt.nullspace()                       # columns of X span ker Mt (over Q), nul = dimension
assert nul == U - rk
B = [[int(X[i, j]) for i in range(U)] for j in range(nul)]
# make primitive integer rows
Bp = []
for v in B:
    g = 0
    for t in v: g = math.gcd(g, t)
    Bp.append([t // g for t in v] if g else v)
Bm = flint.fmpz_mat(Bp)
print(f"nullspace basis: max |entry| = {max(abs(t) for v in Bp for t in v)}  ({time.time()-t0:.0f}s)", flush=True)
# saturate: the Z-span of a Q-basis need not be V cap Z^A; LLL of the basis of the saturated lattice: use hnf-free approach:
# lattice V cap Z^A = kernel lattice; flint nullspace returns rational-reduced basis; we LLL it and then check saturation
# by the gcd of maximal minors being 1 is too costly; we report heights for the LLL basis of the lattice spanned (an
# upper bound for the heights needed, since that lattice is contained in V cap Z^A).
Lr = Bm.lll()
red = np.array([[int(Lr[i, j]) for j in range(U)] for i in range(nul)], dtype=object)
heights = sorted(set(int(max(abs(int(t)) for t in row)) for row in red))
Hred = heights[-1]
l1 = [sum(abs(int(t)) for t in row) for row in red]
print(f"LLL basis: H_red = max|c|_inf = {Hred}; height profile (H: #vectors) = "
      f"{[(h, sum(1 for row in red if max(abs(int(t)) for t in row) <= h)) for h in heights]}; max |c|_1 = {max(l1)}  ({time.time()-t0:.0f}s)", flush=True)
Af = A.astype(np.float64); redf = red.astype(np.float64)
for th in (0.0, np.pi / (2 * q)):
    u = np.cos(th + 2 * np.pi * np.arange(q) / q); ui = np.abs(u).max()
    g = ((Af @ u) / ui) ** m                                 # hat u
    phis = redf @ g
    kap_single = max(abs(phis[i]) / (l1[i] * q ** m) for i in range(nul))
    out = []
    for h in heights:
        sel = [i for i in range(nul) if max(abs(int(t)) for t in red[i]) <= h]
        Qh, _ = np.linalg.qr(redf[sel].T)
        pr = Qh.T @ g
        psi_h = pr @ pr / U
        out.append((h, math.sqrt(psi_h) / q ** m))
    print(f"theta={th:.4f}: best single LLL vector kappa_c(hat u) = {kap_single:.6e};  lower bounds kappa_H(hat u) >= sqrt(Psi_span)/d^m:"
          f" {[(h, float(f'{v:.6e}')) for h, v in out]}", flush=True)
print(f"done ({time.time()-t0:.0f}s)")
