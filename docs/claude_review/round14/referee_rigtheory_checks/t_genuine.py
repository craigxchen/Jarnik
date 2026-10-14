import sys, random
from glib import *
ps = [5, 13, 17, 29]; pis = [two_sq(p) for p in ps]
V = VK(4, 5); bad = 0; cnt = 0
for q in primerange(3, 400):
    if q in ps: continue
    m, L = genuine(ps, q); n = m//4
    assert [x % 2 for x in L] == nu(ps, q), q
    for v in V:
        z = Zv(pis, v); t = sum(a*b for a, b in zip(v, L)) % m
        for s in range(4):
            cnt += 1
            if (t == s*n) != (Ds(z, s) % q == 0): bad += 1
print("checks", cnt, "mismatches", bad)
