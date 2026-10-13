"""Independent randomized search for the K=4 'flexible top' (ps={5,13,17,29}): genuine labels (all levels) at forced
q <= 103, collision-free above 199; at each of the 18 forced moduli in [107,199] choose a parity-respecting level-1
vector OUTSIDE the genuine class that fits the remaining per-unit budget; then choose level-2 lifts avoiding every
level-2 collision (so deeper levels are clean). Then test system-level coherence on the top Q_7 moduli
{127,151,167,191,199} against all pairs (zeta_a, zeta_b), Norm <= 2000, for two column pairs."""
import sys, math, random, itertools
import numpy as np
from rlib import *
ps = [5, 13, 17, 29]; r = 4; K = 4
F = [3, 7, 11, 19, 23, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97, 101, 103, 107, 109, 113, 127,
     131, 137, 139, 149, 151, 157, 163, 167, 173, 179, 181, 191, 193, 199]
low = [q for q in F if q <= 103]; top = [q for q in F if q > 103]
V = VK_mod_sign(r, K); Varr = np.array(V, dtype=np.int64)
pis = gaussian_primes(ps); Zs = [Zv(pis, v) for v in V]
def vq(x, q):
    x = abs(x); e = 0
    while x and x % q == 0: x //= q; e += 1
    return e
# capacity in log units, exact ints: need (prod)^2 <= num
num = np.zeros((len(V), 4), dtype=object); used = np.ones((len(V), 4), dtype=object)
for i, z in enumerate(Zs):
    nz = z[0] ** 2 + z[1] ** 2
    for s in range(4):
        D = Ds(z, s); pr = 1
        for q in low: pr *= q ** vq(D, q)
        num[i, s] = nz * (2 if s % 2 else 1); used[i, s] = pr
seed = int(sys.argv[1]) if len(sys.argv) > 1 else 1
rng = np.random.default_rng(seed)
chosen = {}
for q in top:
    m = q - chi4(q); n = m // 4; nv = nu(ps, q)
    T = Torus(q); io = T.iota(); L = [T.ratio_log(pi) for pi in pis]
    if io != m // 4: L = [(-x) % m for x in L]
    units = [u for u in range(1, n) if math.gcd(u, n) == 1]
    cls = set(tuple((u * x) % n for x in L) for u in units)
    # allowed[i,s]: one more collision at q fits
    allowed = np.array([[ (used[i, s] * q) ** 2 <= num[i, s] for s in range(4)] for i in range(len(V))])
    found = None
    for batch in range(200):
        g = rng.integers(0, m // 2, size=(100000, r)) * 2 + np.array(nv)[None, :]
        t = (g @ Varr.T) % m
        coll = (t % n == 0); sv = (t // n) % 4
        bad = coll & ~allowed[np.arange(len(V))[None, :], sv]
        ok = np.nonzero(~bad.any(axis=1))[0]
        for idx in ok:
            if tuple(g[idx] % n) not in cls:
                found = g[idx]; break
        if found is not None: break
    if found is None:
        print(f"q={q}: no non-genuine-class fitting vector found in 2e7 samples"); sys.exit(0)
    l = [int(x) for x in found]
    tt = [(sum(a * b for a, b in zip(v, l))) % m for v in V]
    colliding = [(i, (tt[i] // n) % 4) for i in range(len(V)) if tt[i] % n == 0]
    for (i, s) in colliding: used[i, s] *= q
    # level-2 lift: t in (Z/q)^r with q not dividing d_v + 4<v,t> for each colliding v (d_v = <v,l>/n as integer)
    ok2 = False
    for _ in range(2000):
        tv = [int(x) for x in rng.integers(0, q, size=r)]
        if all((sum(a * b for a, b in zip(V[i], l)) // n + 4 * sum(a * b for a, b in zip(V[i], tv))) % q != 0
               for (i, s) in colliding):
            ok2 = True; break
    chosen[q] = (l, len(colliding), ok2)
    print(f"q={q} n={n}: non-genuine-class level-1 vector {l}, {len(colliding)} level-1 collisions, level-2 clean lift found: {ok2}", flush=True)
# final budget check (level 1 at top + all levels genuine low)
viol = sum(1 for i in range(len(V)) for s in range(4) if used[i, s] ** 2 > num[i, s])
print("budget violations:", viol)
# coherence on top Q7
H = 2000
zetas = [(a, b) for a in range(1, 45) for b in range(0, 45) if a * a + b * b <= H]
topQ7 = [q for q in top if q % 8 == 7]
th_forced_Q7 = sum(math.log(q) for q in F if q % 8 == 7)
print("top Q7", topQ7, "share of forced Q7 mass %.3f" % (sum(math.log(q) for q in topQ7) / th_forced_Q7))
for pr in [(0, 1), (2, 3), (0, 2), (1, 3)]:
    Z = len(zetas); acc = np.zeros((Z, Z)); allc = np.ones((Z, Z), bool)
    for q in topQ7:
        n = (q + 1) // 4; T = Torus(q); l = chosen[q][0]; c = [x % n for x in l]
        units = np.array([u for u in range(1, n) if math.gcd(u, n) == 1])
        Lz = np.array([T.ratio_log(z) % n if (z[0] ** 2 + z[1] ** 2) % q else -1 for z in zetas])
        P = (Lz[:, None] * units[None, :]) % n
        A = (P == c[pr[0]]).astype(np.float32); A[Lz < 0] = 1
        B = (P == c[pr[1]]).astype(np.float32); B[Lz < 0] = 1
        coh = (A @ B.T) > 0
        allc &= coh; acc += coh * math.log(q)
    print(f"columns {pr}: zeta pairs coherent on ALL top Q7 moduli: {int(allc.sum())}; max coherent share of top Q7 mass {acc.max()/sum(math.log(q) for q in topQ7):.3f}")
