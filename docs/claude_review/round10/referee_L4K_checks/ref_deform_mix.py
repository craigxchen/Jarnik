"""Independent one-modulus deformation test (level-1 relaxation at q0), K=4, r=4, ps={5,13,17,29}.
System: genuine labels (all levels) at every forced q != q0, collision-free above the forced range (zero demand),
arbitrary parity-respecting l at q0. l survives iff every level-1 collision (v,s) at q0 fits:
  (q0 * prod_{q in F minus q0} q^{v_q(D_s(Z_v))})^2 <= Norm(Z_v) * (2 if s odd).
Reports #survivors and #survivors outside the genuine class c = u L mod n."""
import sys, math, itertools
import numpy as np
from rlib import *
ps = [5, 13, 17, 29]; r = 4; K = int(sys.argv[1]) if len(sys.argv) > 1 else 4
F = {2: [3, 7, 11, 19, 23, 31, 41],
     3: [3, 7, 11, 19, 23, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 83, 103],
     4: [3, 7, 11, 19, 23, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97, 101, 103, 107, 109, 113, 127,
         131, 137, 139, 149, 151, 157, 163, 167, 173, 179, 181, 191, 193, 199],
     5: [3, 7, 11, 19, 23, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97, 101, 103, 107, 109, 113, 127,
         131, 137, 139, 149, 151, 157, 163, 167, 173, 179, 181, 191, 193, 197, 199, 211, 223, 227, 229, 233, 239, 241,
         251, 257, 263, 269, 271, 277, 281, 283, 293, 307, 311, 313, 317, 331, 337, 347, 349]}[int(sys.argv[2])]
q0s = [int(x) for x in sys.argv[3:]]
V = VK_mod_sign(r, K); Varr = np.array(V, dtype=np.int64)
pis = gaussian_primes(ps)
Zs = [Zv(pis, v) for v in V]
def vq(x, q):
    x = abs(x); e = 0
    while x and x % q == 0: x //= q; e += 1
    return e
for q0 in q0s:
    m = q0 - chi4(q0); n = m // 4; nv = nu(ps, q0)
    # remaining capacity check: allowed[i,s] True iff a collision of (v_i,s) at q0 fits
    allowed = np.zeros((len(V), 4), bool)
    for i, z in enumerate(Zs):
        nz = z[0] ** 2 + z[1] ** 2
        for s in range(4):
            D = Ds(z, s)
            pr = q0
            for q in F:
                if q != q0: pr *= q ** vq(D, q)
            allowed[i, s] = (pr * pr <= nz * (2 if s % 2 else 1))
    T = Torus(q0); io = T.iota()
    L = [T.ratio_log(pi) for pi in pis]
    if io != m // 4: L = [(-x) % m for x in L]
    units = [u for u in range(1, max(n, 2)) if math.gcd(u, n) == 1] if n > 1 else [0]
    cls = set(tuple((u * x) % n for x in L) for u in units) if n > 1 else None
    choices = [np.array([x for x in range(m) if x % 2 == nv[j]], dtype=np.int64) for j in range(r)]
    surv = 0; surv_out = 0; total = 0
    A0 = choices[0]
    for x1 in choices[1]:
        grid = np.array(np.meshgrid(A0, [x1], choices[2], choices[3], indexing='ij')).reshape(4, -1).T
        t = (grid @ Varr.T) % m
        coll = (t % n == 0); s = (t // n) % 4
        bad = coll & ~allowed[np.arange(len(V))[None, :], s]
        ok = ~bad.any(axis=1)
        total += len(grid)
        if ok.any():
            g = grid[ok]
            surv += len(g)
            if cls is not None:
                c = g % n
                inc = np.array([tuple(row) in cls for row in c])
                surv_out += int((~inc).sum())
    print(f"monomial norm<={K}, genuine on forced set of norm {sys.argv[2]}: q0={q0} n={n} parity vectors={total} level-1 survivors={surv} outside genuine class={surv_out}", flush=True)
