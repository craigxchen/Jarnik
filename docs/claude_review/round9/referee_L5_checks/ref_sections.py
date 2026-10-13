#!/usr/bin/env python3
"""Referee check of the evidence in L5.md Prop 5.2 ('generic profiles seem to have dense interfaces').

sigma(A) = min over INTEGRAL sections s (A s = I_M, s in Z^(r x M)) of max_j sum_x |s_jx|.
The write-up only reports sections A_S^(-1) from square unimodular minors (medians 5,5,25,19,407,532).
Here sigma(A) is computed by a MILP (scipy/HiGHS) over all integral sections, for random 0/1 profiles
with distinct rows and no constant column; 'opt' = proven optimum, otherwise best found and LP bound.
"""
import random, math, sys
import numpy as np
from scipy.optimize import milp, LinearConstraint, Bounds
random.seed(17)

def sigma_milp(A, tlim):
    M, r = len(A), len(A[0])
    # variables: sp[j,x], sm[j,x] (j<r, x<M), then t ; index (j*M + x)
    nv = 2 * r * M + 1
    c = np.zeros(nv); c[-1] = 1.0
    rows = []; lo = []; hi = []
    for y in range(M):            # (A s)_{y x} = delta_{yx}
        for x in range(M):
            row = np.zeros(nv)
            for j in range(r):
                if A[y][j]:
                    row[j * M + x] = A[y][j]; row[r * M + j * M + x] = -A[y][j]
            rows.append(row); lo.append(1.0 if x == y else 0.0); hi.append(1.0 if x == y else 0.0)
    for j in range(r):            # sum_x sp + sm <= t
        row = np.zeros(nv)
        for x in range(M):
            row[j * M + x] = 1; row[r * M + j * M + x] = 1
        row[-1] = -1
        rows.append(row); lo.append(-np.inf); hi.append(0.0)
    cons = LinearConstraint(np.array(rows), lo, hi)
    integ = np.ones(nv); integ[-1] = 0
    ub = np.full(nv, 60.0); ub[-1] = np.inf
    res = milp(c, constraints=cons, integrality=integ, bounds=Bounds(np.zeros(nv), ub),
               options={'time_limit': tlim, 'disp': False})
    lp = milp(c, constraints=cons, integrality=np.zeros(nv), bounds=Bounds(np.zeros(nv), ub))
    best = res.fun if res.x is not None else None
    if res.x is not None:
        s = np.round(res.x[:r * M]).reshape(r, M) - np.round(res.x[r * M:2 * r * M]).reshape(r, M)
        assert np.array_equal(np.array(A) @ s, np.eye(M)), 'section check'
        best = int(np.abs(s).sum(axis=1).max())
    return best, res.status == 0, lp.fun

print(' M   r   per-instance sigma(A) [opt/* = time-limited best] ; LP lower bound')
for (M, r, k, tl) in [(8, 16, 6, 30), (12, 24, 5, 60), (16, 32, 5, 90), (16, 64, 4, 90)]:
    out = []
    got = 0
    while got < k:
        A = [[random.randint(0, 1) for _ in range(r)] for _ in range(M)]
        if len(set(map(tuple, A))) < M: continue
        if any(len(set(A[i][j] for i in range(M))) == 1 for j in range(r)): continue
        if np.linalg.matrix_rank(np.array(A)) < M: continue
        got += 1
        best, opt, lpv = sigma_milp(A, tl)
        out.append('%s%s(LP %.2f)' % (best, '' if opt else '*', lpv))
    print('%2d %3d   %s' % (M, r, '  '.join(out)))
    sys.stdout.flush()
