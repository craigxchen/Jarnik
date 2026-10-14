"""Exact check: genuine labels collide at q (level 1) with unit s  iff  q | D_s(Z_v)."""
from rlib14 import *
bad = 0; tot = 0
for ps in ([5, 13, 17, 29], [37, 41, 53, 61], [5, 13, 17, 29, 37]):
    pis = gaussian_primes(ps)
    V = VK(len(ps), 3)
    for q in primes_upto(150):
        if q == 2 or any(p % q == 0 for p in ps): continue
        m, L = genuine_labels(ps, q)
        n = m // 4
        nv = nu(ps, q)
        assert all(L[j] % 2 == nv[j] for j in range(len(ps))), (ps, q)
        for v in V:
            t = sum(a * b for a, b in zip(v, L)) % m
            z = Zv(pis, v)
            for s in range(4):
                lhs = (t == s * n)
                rhs = (Ds(z, s) % q == 0)
                tot += 1; bad += (lhs != rhs)
print("checks", tot, "mismatches", bad)
