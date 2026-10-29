#!/usr/bin/env python3
"""Cor 3.4 test (B''): random (mostly non-nef) torus-invariant D with D.beta = 8t; dim A_t(P_D cap M) <= 1 ?"""
import random, sys, time
sys.path.insert(0, '.')
from cor34 import PD_points, RAYS
from thm33 import dimA
rnd = random.Random(7)
hist = {}; t0 = time.time(); n = 0
for trial in range(20000):
    d = {A: 0 for A in RAYS}
    for _ in range(rnd.choice([8, 16, 24])):
        d[rnd.choice(RAYS)] += 1
    # allow a few negative coefficients (non-effective D, still possibly effective class)
    for _ in range(rnd.randint(0, 3)):
        A = rnd.choice(RAYS); d[A] -= 1; B = rnd.choice(RAYS); d[B] += 1
    Db = sum(d.values()); t = Db // 8
    S = PD_points(d)
    if not S or len(S) > 400: continue
    dm = dimA(S, t)
    hist[dm] = hist.get(dm, 0) + 1; n += 1
print(f"(B'') {n} random classes with nonempty P_D and D.beta = 8t: histogram of h^0(pi*D - tE): {hist}  [{time.time()-t0:.1f}s]")
