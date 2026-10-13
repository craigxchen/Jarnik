#!/usr/bin/env python3
"""Evidence: independent uniform residues at odd prime powers <= X = 3M (random lifts), demand D,
lambda_max(D)/(M log M), max_x D_x/(M log M), D_max/log M. Compare with p3lower Sec. 2 table."""
import math, numpy as np
from sympy import primerange
rng = np.random.default_rng(31)
def mq(q): return q - 1 if q % 4 == 1 else q + 1
for M in [64, 256, 1024]:
    X = 3 * M
    D = np.zeros((M, M))
    for q in primerange(3, X + 1):
        Q = q; m = mq(q); kap = rng.integers(0, m, size=M)
        while Q <= X:
            eq = kap[:, None] == kap[None, :]
            D += 2 * math.log(q) * eq
            Q *= q
            if Q <= X:
                kap = kap + m * rng.integers(0, q, size=M); m *= q
    np.fill_diagonal(D, 0)
    L = math.log(M)
    print(f"M={M}: lambda_max/(M log M)={np.linalg.eigvalsh(D)[-1]/(M*L):.3f}, max_x D_x/(M log M)={D.sum(1).max()/(M*L):.3f}, D_max/log M={D.max()/L:.2f}")
