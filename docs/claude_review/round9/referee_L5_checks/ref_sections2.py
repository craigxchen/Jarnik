#!/usr/bin/env python3
"""Upper bounds on sigma(A) at M = 12, 16 (supplement to ref_sections.py).
For each row index x solve the small MILP  min ||s_x||_1  s.t.  A s_x = e_x, s_x integral  (scipy/HiGHS),
assemble s = (s_x)_x (an integral section, checked exactly) and report max_j sum_x |s_jx| >= sigma(A).
Same random-profile distribution as ref_sections.py (0/1, distinct rows, no constant column, rank M)."""
import random, sys
import numpy as np
from scipy.optimize import milp, LinearConstraint, Bounds
random.seed(23)

def col_section(A, x, tlim=20):
    M, r = len(A), len(A[0])
    Aeq = np.hstack([np.array(A, float), -np.array(A, float)])
    b = np.zeros(M); b[x] = 1
    res = milp(np.ones(2 * r), constraints=LinearConstraint(Aeq, b, b), integrality=np.ones(2 * r),
               bounds=Bounds(np.zeros(2 * r), np.full(2 * r, 50.0)), options={'time_limit': tlim})
    if res.x is None: return None
    v = np.round(res.x[:r]) - np.round(res.x[r:])
    return v

print(' M   r   upper bounds on sigma(A) (row-l1 of assembled column-wise l1-minimal integral section)')
for (M, r, k) in [(12, 24, 8), (16, 32, 8), (16, 64, 4), (24, 48, 3)]:
    out = []; got = 0
    while got < k:
        A = [[random.randint(0, 1) for _ in range(r)] for _ in range(M)]
        if len(set(map(tuple, A))) < M: continue
        if any(len(set(A[i][j] for i in range(M))) == 1 for j in range(r)): continue
        if np.linalg.matrix_rank(np.array(A)) < M: continue
        got += 1
        cols = [col_section(A, x) for x in range(M)]
        if any(c is None for c in cols):
            out.append('none'); continue
        s = np.array(cols).T
        assert np.array_equal((np.array(A) @ s).astype(int), np.eye(M, dtype=int))
        out.append(str(int(np.abs(s).sum(axis=1).max())))
    print('%2d %3d   %s' % (M, r, ' '.join(out))); sys.stdout.flush()
