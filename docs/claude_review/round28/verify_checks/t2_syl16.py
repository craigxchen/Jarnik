"""Theory §2.3 check: Sylvester n=16 (d=15), T=1, L=2,3,4: h = max sigma on A, number of top points,
(R) = no top point has a legal pair move k + e_x - e_y in A; also max sigma on A_{e_0-e_1}.  Exact enumeration of 3^15."""
import sys, itertools
import numpy as np
sys.path.insert(0, '/tmp/claude-0/-home-user-Jarnik/ba4929a4-8665-5790-985a-3f09130a14b9/scratchpad/round28/verify_checks')
from vlib import *

n, pats = syl(4); d = 15
# functional matrix: kfull(S) as a functional on k' (base = 15)
F = np.array([[(1 if x in S else 0) - (1 if d in S else 0) for x in range(d)] for S in pats], dtype=np.int16)
suf = np.array(list(itertools.product([-1, 0, 1], repeat=10)), dtype=np.int16)
for L in (2, 3, 4):
    best = -99; tops = []
    best01 = -99
    for pre in itertools.product([-1, 0, 1], repeat=5):
        P = np.concatenate([np.tile(np.array(pre, dtype=np.int16), (len(suf), 1)), suf], axis=1)
        ok = (np.abs(P @ F.T) <= L).all(axis=1)
        Q = P[ok]
        if len(Q) == 0: continue
        s = Q.sum(axis=1)
        m = s.max()
        if m > best: best = m; tops = [Q[s == m]]
        elif m == best: tops.append(Q[s == m])
        # A_{e0-e1}: k in A with k + e0 - e1 in A
        Q2 = Q.copy(); Q2[:, 0] += 1; Q2[:, 1] -= 1
        ok2 = (np.abs(Q2) <= 1).all(axis=1) & (np.abs(Q2 @ F.T) <= L).all(axis=1)
        if ok2.any(): best01 = max(best01, Q[ok2].sum(axis=1).max())
    T = np.concatenate(tops)
    moves = 0
    for x in range(d):
        for y in range(d):
            if x == y: continue
            M = T.copy(); M[:, x] += 1; M[:, y] -= 1
            okm = (np.abs(M) <= 1).all(axis=1) & (np.abs(M @ F.T) <= L).all(axis=1)
            moves += int(okm.sum())
    print('syl16 L', L, 'h', int(best), 'averaging bound floor(2L(1-1/16))', (2 * L * 15) // 16, 'top points', len(T),
          'legal pair moves', moves, '(R)', moves == 0, 'max sigma on A_{e0-e1}', int(best01), flush=True)
