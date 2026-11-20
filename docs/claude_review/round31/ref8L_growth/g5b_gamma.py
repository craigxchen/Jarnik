#!/usr/bin/env python3
"""Search squarefree N (products of split primes) for arcs of angular width C N^{-1/4} (C=2,3) whose marks have
Gaussian gcd 1 and contain a pair of primitive height H close to N^{1/2}. Report min log H/log N."""
import math, random, itertools, sys
sys.path.insert(0, '.')
from g1_arcs import points, clusters, ggcd, gnorm, gdivmod_exact
from sympy import primerange
random.seed(5)
split = [p for p in primerange(5, 3000) if p % 4 == 1]
best = {2.0: (9, None), 3.0: (9, None)}
for trial in range(3000):
    k = random.randint(3, 6)
    ps = random.sample(split, k)
    fac = {p: 1 for p in ps}
    N = math.prod(ps)
    _, pts = points(fac)
    for C in (2.0, 3.0):
        for cl in clusters(pts, C*N**-0.25):
            if len(cl) < 3: continue
            zs = [t[0] for t in cl]
            g = zs[0]
            for z in zs[1:]: g = ggcd(g, z)
            if gnorm(g) != 1: continue
            for z, w in itertools.combinations(zs, 2):
                H = N // gnorm(ggcd(z, w))
                r = math.log(H)/math.log(N)
                if r < best[C][0]: best[C] = (r, (N, len(cl)))
print(best)
