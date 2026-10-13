# Joint search on ALL forced moduli (K=4 tight set): can a good system be non-genuine-class on the core?
import sys, math, random, time
import numpy as np
import forced_coh
from forced_coh import *
K, r, tag, qcore, tries, cap = int(sys.argv[1]), int(sys.argv[2]), sys.argv[3], int(sys.argv[4]), int(sys.argv[5]), int(sys.argv[6])
P, F = load_set(K, r, tag); V = monos(r, K)
B = np.array([[budget(P, v, s) for s in range(4)] for v in V])
mods = [Mod(q, P, K, V) for q in F]
def isgen(md, ell):
    n = md.n; Lc = [x % n for x in md.L]
    return any(all((u * Lc[j] - ell[j]) % n == 0 for j in range(len(ell))) for u in range(1, n + 1) if math.gcd(u, n) == 1)
found = 0; t0 = time.time(); depths = []
for seed in range(tries):
    rng = random.Random(1000 + seed)
    g = gen_good2(mods, B, rng, per_q=6, node_cap=cap)
    depths.append(gen_good2.depth)
    if g is None: continue
    found += 1; sysd = g[0]
    ng = [q for q in F if not isgen(forced_coh_mod := [md for md in mods if md.q == q][0], sysd[q][0])]
    print(f"  seed {seed}: good system on all {len(F)} forced moduli; non-genuine-class at core q<={qcore}: {[q for q in ng if q <= qcore]}; at top: {len([q for q in ng if q > qcore])}", flush=True)
print(f"K={K} r={r} {tag}: {found}/{tries} searches found a good system on all forced moduli; max depth reached per search {depths} of {len(F)} [{time.time()-t0:.0f}s]")
