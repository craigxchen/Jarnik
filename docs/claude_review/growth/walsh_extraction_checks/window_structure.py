"""Exact parallelograms (F_2 additive quadruples on the varying primes) inside windows of actual circles:
extremal windows versus typical windows with the same point count, versus random subsets.
Classifies each parallelogram by its matching weights (12|34, 13|24, 14|23) and detects integer rectangles
(z1 z2 = unit * z3 z4 exactly) and affine 3-cubes (pure order-8 Walsh cores)."""
import numpy as np, math, random, sys, os, itertools
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from gauss_common import *
from cluster_search import angles_for, bits_of, SP

Q = math.pi / 2

def window_sets(P, M, nrand=40):
    pis, ang = angles_for(P)
    a = np.mod(ang, Q); order = np.argsort(a); srt = a[order]; n = len(srt)
    ext = np.concatenate([srt, srt + Q])
    widths = ext[M - 1:M - 1 + n] - ext[:n]
    i = int(np.argmin(widths))
    best = [int(order[(i + t) % n]) for t in range(M)]
    rnd = []
    for _ in range(nrand):
        j = random.randrange(n)
        rnd.append([int(order[(j + t) % n]) for t in range(M)])
    return best, rnd, widths[i]

def analyse(P, idxs):
    k = len(P)
    S = np.array([[(-1) ** b for b in bits_of(ix, k)] for ix in idxs])
    w = np.array([math.log(p) for p in P])
    var = [j for j in range(k) if len(set(S[:, j])) > 1]
    Sv = S[:, var]; wv = w[var]; W = wv.sum()
    M = len(idxs)
    npar = 0; nrect = 0; mins12 = []
    for Qd in itertools.combinations(range(M), 4):
        r = [Sv[q] for q in Qd]
        if np.all(r[0] * r[1] * r[2] * r[3] > 0):
            npar += 1
            # matching weights
            m = {}
            for (A, B) in [((0, 1), (2, 3)), ((0, 2), (1, 3)), ((0, 3), (1, 2))]:
                cols = (r[A[0]] == r[A[1]]) & (r[B[0]] == r[B[1]]) & (r[A[0]] != r[B[0]])
                m[(A, B)] = wv[cols].sum() / W
            vals = sorted(m.values())
            mins12.append(vals[0])
            if vals[0] < 1e-12: nrect += 1   # some matching absent => integer rectangle (multiplicative quadruple)
    # affine 3-cubes
    ncube = 0
    if M <= 16:
        X = [tuple(int(v < 0) for v in row) for row in Sv]
        Xs = set(X)
        for sub in itertools.combinations(range(M), 8):
            base = X[sub[0]]
            D = set(tuple(a ^ b for a, b in zip(X[s], base)) for s in sub)
            if all(tuple(a ^ b for a, b in zip(u, v)) in D for u in D for v in D):
                ncube += 1
    return npar, nrect, ncube, (min(mins12) if mins12 else None)

random.seed(3)
tot = {}
for trial in range(60):
    k = random.randint(12, 17)
    P = sorted(random.sample(SP[:k + 20], k))
    if sum(math.log10(p) for p in P) < 18: continue
    for M in (8, 12, 16):
        best, rnd, wd = window_sets(P, M, nrand=15)
        b = analyse(P, best)
        rr = [analyse(P, r) for r in rnd]
        # random subsets of the cube
        rs = [analyse(P, random.sample(range(2 ** k), M)) for _ in range(15)]
        t = tot.setdefault(M, [0, 0, 0, 0, 0, 0, 0, 0, 0, 0])
        t[0] += 1; t[1] += b[0]; t[2] += b[1]; t[3] += b[2]
        t[4] += np.mean([x[0] for x in rr]); t[5] += np.mean([x[1] for x in rr]); t[6] += np.mean([x[2] for x in rr])
        t[7] += np.mean([x[0] for x in rs]); t[8] += np.mean([x[1] for x in rs]); t[9] += np.mean([x[2] for x in rs])
for M, t in sorted(tot.items()):
    n = t[0]
    print("M=%2d circles=%2d | extremal window: par=%.2f rect=%.2f cube8=%.2f | typical window: par=%.2f rect=%.2f cube8=%.3f | random subset: par=%.3f rect=%.3f cube8=%.3f"
          % (M, n, t[1] / n, t[2] / n, t[3] / n, t[4] / n, t[5] / n, t[6] / n, t[7] / n, t[8] / n, t[9] / n))
