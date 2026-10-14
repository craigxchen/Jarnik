"""Check: genuine labels collide at q^a with unit s iff q^a | D_s(Z_v); and parity L_j = nu_j mod 2."""
import sys
from rlib import *
sets = [[5,13,17,29],[173,181,193,197],[37,41,53,61,73]]
nchk = 0; bad = 0; parbad = 0
for ps in sets:
    pis = gaussian_primes(ps)
    N = 1
    for p in ps: N *= p
    for q in primes_upto(110):
        if q == 2 or N % q == 0: continue
        for a in (1, 2):
            if a == 2 and q > 60: continue
            T = Torus(q, a)
            L = [T.ratio_log(pi) for pi in pis]
            io = T.iota(); m = T.order
            if a == 1:
                nv = nu(ps, q)
                for j in range(len(ps)):
                    if L[j] % 2 != nv[j]: parbad += 1
            for v in VK(len(ps), 3):
                z = Zv(pis, v)
                t = sum(x * l for x, l in zip(v, L)) % m
                for s in range(4):
                    lab = (t == (s * io) % m)
                    ar = (Ds(z, s) % (q ** a) == 0)
                    nchk += 1
                    if lab != ar: bad += 1
print("checks", nchk, "mismatches", bad, "parity mismatches", parbad)
