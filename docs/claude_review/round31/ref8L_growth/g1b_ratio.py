#!/usr/bin/env python3
"""Targeted test of the ratio multiplicity bound (VIII.5): N = n (a^2+1), u=(a+i)/(a-i).
For each arc (window starting at each mark) of angular width C N^{-1/4} count ordered pairs per ratio;
report the minimal C at which a multiplicity r is attained and compare with r <= 1 + C^2."""
import math, itertools, sys
from sympy import factorint
sys.path.insert(0, '.')
from g1_arcs import points, arg, reduce_ratio, normalize_unit, gnorm
worst = 0.0
for n in (25, 65, 325, 1105, 5525):
    for a in range(2, 400, 2):
        fac = factorint(n*(a*a+1))
        if any(p % 4 != 1 for p in fac): continue
        N = n*(a*a+1)
        ps, pts = points(fac)
        pts = sorted(pts, key=lambda t: arg(t[0]))
        angs = [arg(t[0]) for t in pts]
        L = len(pts)
        # for each pair of marks, required C for both to be in one arc is their angular distance * N^{1/4}.
        # For each ratio, multiplicity among pairs inside a window: brute force over windows.
        for i in range(L):
            # window starting at i, grow
            ratio = {}
            for k in range(1, min(L, 60)):
                w = (angs[(i+k) % L] - angs[i]) % (2*math.pi)
                C = w * N**0.25
                if C > 6: break
                cl = [pts[(i+j) % L][0] for j in range(k+1)]
                newz = cl[-1]
                for z in cl[:-1]:
                    for (x, y) in ((z, newz), (newz, z)):
                        key = normalize_unit(*reduce_ratio(x, y))
                        ratio[key] = ratio.get(key, 0) + 1
                        r = ratio[key]
                        if r > 1:
                            worst = max(worst, (r - 1) / (C*C))
                            if r - 1 > C*C: print('VIOLATION', n, a, r, C)
print('max over tests of (multiplicity-1)/C^2 =', worst)
