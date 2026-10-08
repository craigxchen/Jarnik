"""Like cluster_search.py but only on circles with log10 N >= LMIN (arc a tiny fraction of the circle).
Also records the Poisson expectation lambda = (#points) * C / (2 pi sqrt R) of a window of the same size."""
import numpy as np, math, random, json, sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from gauss_common import *
from cluster_search import best_windows, cluster_points, verify, SP

seed = int(sys.argv[1]); LMIN = float(sys.argv[2]); KMAX = int(sys.argv[3]); NTR = int(sys.argv[4])
random.seed(seed)
Ms = [4, 5, 6, 7, 8, 10, 12, 16]
best = {M: None for M in Ms}
cnt = 0
while cnt < NTR:
    k = random.randint(12, KMAX)
    pool = SP[:random.randint(k, min(60, k + 30))]
    P = sorted(random.sample(pool, k))
    if sum(math.log10(p) for p in P) < LMIN: continue
    cnt += 1
    pis, ang, res = best_windows(P, Ms)
    for M, (C, idxs, width) in res.items():
        if best[M] is None or C < best[M][0]:
            best[M] = (C, P, idxs)
out = {}
for M in Ms:
    C, P, idxs = best[M]
    pis, ang, _ = best_windows(P, [M])
    cl = cluster_points(P, pis, ang, idxs)
    ok, zs = verify(P, pis, cl, C)
    k = len(P); logN = sum(math.log(p) for p in P)
    lam = 4 * 2 ** k * C / (2 * math.pi * math.exp(logN / 4))
    out[M] = {"C": C, "P": P, "signs": [s for (s, u, th) in cl], "units": [u for (s, u, th) in cl],
              "points": [list(z) for z in zs], "verified": ok, "lambda": lam}
    print("M=%3d  C=%9.3f  k=%2d  log10N=%5.1f  lambda=%.3g  verified=%s" % (M, C, k, logN / math.log(10), lam, ok))
json.dump(out, open("large_clusters_s%d_L%d.json" % (seed, int(LMIN)), "w"))
