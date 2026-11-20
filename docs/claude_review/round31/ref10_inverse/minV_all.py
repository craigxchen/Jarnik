#!/usr/bin/env python3
"""Companion to minV.py: (a) q = 31, 43 owner minimum; (b) minimum over V of the weight of ALL rows
(owners and non-owner labels), i.e. sum_{j in I} sum_l |V_jl|, per column.  Exact verification inside."""
import sys, random
sys.path.insert(0, '.')
from minV import setup, solve_col, parity_lb, run
for q in [int(x) for x in sys.argv[1:]]:
    if q >= 31:
        print(run(q), flush=True)
    I, K, S, a, D, R = setup(q)
    n = len(I); tot = 0; rows = [0]*n
    for c in range(n):
        rhs = [R[r][c] for r in range(len(D))]
        val, x = solve_col(D, rhs, n, list(range(n)))
        tot += val
        for r in range(n): rows[r] += abs(x[r])
    own = sorted(rows[I.index(i)] for i in S); non = sorted(rows[r] for r in range(n) if I[r] not in S)
    print(dict(q=q, allrow_min_total=tot, owner_rows=own, nonowner_rows=non), flush=True)
