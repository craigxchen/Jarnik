"""Evidence for H2: pair demands of the coherent (G = I) design with root repair on the sharp section-4
Paley profile (b copies, canonical twins), at all odd prime powers <= X = 3M.
For each pair x<y: D1 = log m_(xy, unit 1) and Dg = log m^gen_xy (some unit). Prints max and mean / log M.
Points are z_x = prod_j t_j^(a_xj) in T_Q (common factors cancel in ratios)."""
import random, sys, math
import numpy as np
from gauss import *
from check_T3 import profile

def run(q0, b=6, seed=5):
    rng = random.Random(seed)
    A, idx, twins = profile(q0, b); M, r = A.shape
    X = 3*M
    ps = []
    while len(ps) < r:
        p = rng.randrange(M**4, 2*M**4)
        if p % 4 == 1 and is_prime(p): ps.append(p)
    hi = 4*X
    while True:
        zs = gaussian_primes_norm_range(X, hi)
        if len(zs) >= r: zs = zs[:r]; break
        hi *= 2
    D1 = np.zeros((M, M)); Dg = np.zeros((M, M))
    for q in primes_upto(X)[1:]:
        Aq = 1
        while q**(Aq+1) <= X: Aq += 1
        Q = q**Aq; f, _ = sylow2_generator(q, Aq, rng)
        ts = []
        for j, (x, y) in enumerate(zs):
            eps = 0 if legendre(x*x+y*y, q) == legendre(ps[j], q) else 1
            z = (x % Q, y % Q)
            ts.append(mul(mul(z, inv(conj(z, Q), Q), Q), pw(f, eps, Q), Q))
        pts = []
        for xr in range(M):
            z = (1, 0)
            for j in range(r):
                if A[xr, j]: z = mul(z, ts[j], Q)
            pts.append(z)
        lq = math.log(q)
        for x in range(M):
            for y in range(x+1, M):
                u = mul(pts[x], inv(pts[y], Q), Q)
                a1 = ag = 0
                for lev in range(1, Aq+1):
                    Ql = q**lev; ul = (u[0] % Ql, u[1] % Ql)
                    if ul == (1, 0): a1 = lev
                    if pw(ul, 4, Ql) == (1, 0): ag = lev
                D1[x, y] += a1*lq; Dg[x, y] += ag*lq
    iu = np.triu_indices(M, 1)
    lM = math.log(M)
    d1, dg = D1[iu], Dg[iu]
    rowsum = (D1 + D1.T).sum(axis=1)
    lam = np.linalg.eigvalsh(D1 + D1.T).max()
    return M, r, d1.max()/lM, d1.mean()/lM, dg.max()/lM, dg.mean()/lM, lam/(M*lM), rowsum.max()/(M*lM)

if __name__ == "__main__":
    qs = [int(a) for a in sys.argv[1:]] or [43, 59, 83, 107, 131]
    print("  M    r  | unit-1: max  mean | some unit: max mean | lambda_max(D1)/(M logM) | max rowsum/(M logM)")
    for q0 in qs:
        M, r, a, b_, c, d, lam, rs = run(q0)
        print(f"{M:4d} {r:5d} |   {a:6.2f} {b_:6.3f}  |   {c:6.2f} {d:6.3f}     |   {lam:6.3f}              | {rs:6.3f}", flush=True)
