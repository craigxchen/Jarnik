"""Independent joint test of the K=4 core {3..47} (ps={5,13,17,29}): genuine labels (all levels) at forced q in (47,199],
collision-free above; enumerate ALL good level-1 systems on the core (core levels >= 2 not charged) by DFS over
signature tuples; report #tuples, #systems, and whether any is outside the genuine class at some core modulus."""
import sys, math, itertools, time
import numpy as np
from collections import Counter, defaultdict
from rlib import *
ps = [5, 13, 17, 29]; r = 4; K = 4
F = [3, 7, 11, 19, 23, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97, 101, 103, 107, 109, 113, 127,
     131, 137, 139, 149, 151, 157, 163, 167, 173, 179, 181, 191, 193, 199]
core = [q for q in F if q <= 47]; rest = [q for q in F if q > 47]
V = VK_mod_sign(r, K); Varr = np.array(V, dtype=np.int64)
pis = gaussian_primes(ps); Zs = [Zv(pis, v) for v in V]
keys = [(i, s) for i in range(len(V)) for s in range(4)]
def vq(x, q):
    x = abs(x); e = 0
    while x and x % q == 0: x //= q; e += 1
    return e
cap2 = []  # remaining squared capacity: Norm*(2 if odd) / (rest demand)^2  (rational -> keep as (num, den))
num = []; den = []
for (i, s) in keys:
    z = Zs[i]; nz = z[0] ** 2 + z[1] ** 2; D = Ds(z, s)
    pr = 1
    for q in rest: pr *= q ** vq(D, q)
    num.append(nz * (2 if s % 2 else 1)); den.append(pr * pr)
per_q = {}; gclass = {}
for q in core:
    m = q - chi4(q); n = m // 4; nv = nu(ps, q)
    T = Torus(q); io = T.iota(); L = [T.ratio_log(pi) for pi in pis]
    if io != m // 4: L = [(-x) % m for x in L]
    choices = [np.array([x for x in range(m) if x % 2 == nv[j]], dtype=np.int64) for j in range(r)]
    grid = np.array(list(itertools.product(*choices)), dtype=np.int64)
    t = (grid @ Varr.T) % m
    coll = (t % n == 0); sv = (t // n) % 4
    units = [u for u in range(1, max(n, 2)) if math.gcd(u, n) == 1] if n > 1 else [0]
    cls = set(tuple((u * x) % n for x in L) for u in units) if n > 1 else None
    sigs = defaultdict(lambda: [0, 0])  # code -> [count, count outside class]
    for row in range(len(grid)):
        cols = np.nonzero(coll[row])[0]
        code = tuple(sorted(4 * int(i) + int(sv[row, i]) for i in cols))
        out = 0 if cls is None else int(tuple(grid[row] % n) not in cls)
        e = sigs[code]; e[0] += 1; e[1] += out
    per_q[q] = [(code, c, o) for code, (c, o) in sigs.items()]
    print(f"q={q} n={n} vectors={len(grid)} signatures={len(sigs)}", flush=True)
order = sorted(core, key=lambda q: -q)
prod = [1] * len(keys)
res = {"tuples": 0, "systems": 0, "nongen_systems": 0}
t0 = time.time()
def dfs(d, mult, mult_all_gen):
    if d == len(order):
        res["tuples"] += 1; res["systems"] += mult; res["nongen_systems"] += mult - mult_all_gen
        return
    q = order[d]
    for code, c, o in per_q[q]:
        ok = True
        for k in code:
            pq = prod[k] * q
            if pq * pq * den[k] > num[k]: ok = False; break
        if not ok: continue
        for k in code: prod[k] *= q
        dfs(d + 1, mult * c, mult_all_gen * (c - o))
        for k in code: prod[k] //= q
dfs(0, 1, 1)
print(f"core {core}: good signature tuples={res['tuples']} good level-1 systems={res['systems']} "
      f"systems outside the genuine class at >=1 core modulus={res['nongen_systems']}  [{time.time()-t0:.1f}s]")
