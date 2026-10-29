#!/usr/bin/env python3
"""Cor 3.4 test (B'): D = D_{Z1} + D_{Z2} (+ D_{Z3}); P_D can be much larger than Z1+Z2; t = total edges.
Count h^0(pi^*D - tE) = dim A_t(P_D cap M); Cor 3.4 predicts exactly 1 (the product)."""
import random, sys, time
sys.path.insert(0, '.')
from cor34 import zonotope, minimal_D, PD_points, RAYS
from thm33 import dimA
rnd = random.Random(99)
pairs = [(i, j) for i in range(5) for j in range(i + 1, 5)]
fails = 0; tested = 0; bigger = 0; t0 = time.time(); maxS = 0
for trial in range(150):
    parts = rnd.randint(2, 3)
    d = {A: 0 for A in RAYS}; t = 0; Zsum = {(0, 0, 0, 0)}
    for _ in range(parts):
        k = rnd.randint(1, 2)
        Z = zonotope([rnd.choice(pairs) for _ in range(k)])
        # random translate of each factor (changes D by a principal divisor only)
        dz = minimal_D(Z)
        for A in RAYS: d[A] += dz[A]
        t += k
        Zsum = {tuple(a + b for a, b in zip(x, y)) for x in Zsum for y in Z}
    S = PD_points(d)
    if len(S) > 700: continue
    tested += 1; bigger += len(S) > len(Zsum); maxS = max(maxS, len(S))
    dm = dimA(S, t)
    if dm != 1:
        fails += 1; print('unexpected dim', dm, 'parts', parts, 't', t, '|P_D cap M|', len(S))
print(f"(B') {tested} classes D = sum of 2-3 root-product divisors, t = #edges: {fails} with h^0 != 1; "
      f"P_D strictly larger than the product zonotope in {bigger}; max |P_D cap M| = {maxS}; {time.time()-t0:.1f}s")
