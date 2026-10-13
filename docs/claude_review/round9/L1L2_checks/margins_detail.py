"""Detailed margin profile of record clusters and scan clusters (evidence).
Prints sorted margins of the record clusters and, over the scan, K with and without complementary
pairs (pairs with g_xy a unit, mu = W/2)."""
import pickle, sys, itertools
from math import log
import numpy as np
from glib import reduced_chord_data, norm
from known import largest_clusters, normalise

def margins(cl):
    cl, N = normalise(cl)
    W = log(N); out = []
    for x, y in itertools.combinations(range(len(cl)), 2):
        g, beta, s, B = reduced_chord_data(cl[x], cl[y])
        out.append((W / 2 - log(norm(g)), norm(g) == 1, B))
    return W, out

for N, Cm in [(1176852625, 8), (235370525, 5.3), (537606725, 5.9)]:
    W, ms = margins(largest_clusters(N, Cm))
    mus = sorted(m for m, _, _ in ms)
    print(f"N={N} M={int((1+(1+8*len(ms))**0.5)/2)} W/2={W/2:.2f} margins sorted:", [round(m, 2) for m in mus])
    print("   complementary pairs:", sum(c for _, c, _ in ms))

data = pickle.load(open(sys.argv[1], 'rb')) if len(sys.argv) > 1 else []
res = {}
for (N, fac, cl, Cc) in data:
    W, ms = margins(cl)
    mus = np.array([m for m, _, _ in ms]); comp = np.array([c for _, c, _ in ms])
    mean = mus.mean()
    if mean <= 0:
        continue
    K = mus.max() / mean
    K2 = mus[~comp].max() / mean if (~comp).any() else float('nan')
    res.setdefault(len(cl), []).append((K, K2, comp.any()))
for M in sorted(res):
    a = np.array(res[M], dtype=float)
    print(f"scan M={M}: n={len(a)} median K={np.median(a[:,0]):.2f}  median K without complementary pairs={np.nanmedian(a[:,1]):.2f}  frac with complementary pair={a[:,2].mean():.3f}")
