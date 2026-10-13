# EVIDENCE. Supercritical model n=4, r=4, T=4, w=6, m=13 (single_tc.py seed 1 profile): tangent-cone residual tau at the
# sparse directions u = e_x - mean (point x moved alone), for every choice of base point. Prop 5.1 puts the direction in
# the base locus when x != base (2T = 8 < m = 13); with x = base it is not covered.
import random, numpy as np
from engine import System
random.seed(1); np.random.seed(1)
n, r, T, w, m = 4, 4, 4, 6, 13
cols = []
while len(cols) < r:
    a = tuple(random.randint(0, 1) for _ in range(n))
    if 0 < sum(a) < n: cols.append(a)
print("cols", cols)
for b in range(n):
    S = System(list(range(n)), b, cols, T, w, m)
    taus = []
    for x in range(n):
        u = -np.ones(n)/n; u[x] += 1.0
        taus.append(S.tau(u))
    print(f"base={b}: dim V={S.dimV}; tau(e_x - mean) for x=0..3: " + ", ".join(f"{t:.2e}" for t in taus))
