"""Randomised greedy search for lattice zonotopes k0 + sum_i [0,L_i] v_i inside a T = 1 box A (A = lattice points of the
convex polytope cube cap slabs, so containment = all 2^s vertices in the polytope).  Reports the best sum L_i found
(a lower bound for zeta(A)); the found zonotope is re-verified point by point."""
import sys, itertools, random, time
import numpy as np
from flib import *
from model import profile

name = sys.argv[1]; L = int(sys.argv[2]); trials = int(sys.argv[3]) if len(sys.argv) > 3 else 2000
M, pats = profile(name)
d = M - 1
cols = [rel(P) for P in pats]
F = np.zeros((len(cols) + 1, d), dtype=np.int64)
for j, C in enumerate(cols):
    F[j, C] = 1
F[-1, :] = 1  # base point: |sum k'| <= 1
caps = np.array([L] * len(cols) + [1])


def inpoly(K):
    K = np.atleast_2d(K)
    return (np.abs(K) <= 1).all(axis=1) & (np.abs(K @ F.T) <= caps).all(axis=1)


A = box(M, pats, L)
Aset = set(A)
# candidate generators: primitive v with support <= 3, entries in {-1,0,1}, sign-normalised
cands = []
for s in (1, 2, 3, 4):
    for X in itertools.combinations(range(d), s):
        for sg in itertools.product((-2, -1, 1, 2), repeat=s):
            if sg[0] > 0 and (s > 1 or sg[0] == 1):
                v = np.zeros(d, np.int64); v[list(X)] = sg
                cands.append(v)
rng = random.Random(1)
best = (0, None)
t0 = time.time()
for tr in range(trials):
    k0 = np.array(rng.choice(A), dtype=np.int64)
    verts = k0[None, :].copy()
    gens = []
    total = 0
    order = list(range(len(cands)))
    rng.shuffle(order)
    for i in order:
        v = cands[i]
        for Lv in (2, 1):
            nv = np.concatenate([verts, verts + Lv * v[None, :]])
            if inpoly(nv).all():
                verts = nv
                gens.append((tuple(int(x) for x in v), Lv))
                total += Lv
                break
        if len(verts) > 4096:
            break
    if total > best[0]:
        best = (total, (tuple(int(x) for x in k0), gens))
        print('trial %d: sum L = %d  k0=%s gens=%s (%.1fs)' % (tr, total, best[1][0], gens, time.time() - t0), flush=True)
# verify best point by point
k0, gens = best[1]
pts = [np.array(k0)]
for v, Lv in gens:
    pts = [p + t * np.array(v) for p in pts for t in range(Lv + 1)]
assert all(tuple(int(x) for x in p) in Aset for p in pts)
print('RESULT %s L=%d: zonotope with sum L_i = %d verified inside A (|A|=%d)' % (name, L, best[0], len(A)))
