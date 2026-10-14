# Exact check: for genuine labels, v collides at q (level 1) with unit s  <=>  q | D_s(Z_v).
# Collision unit computed in T_q via the discrete log of i (iota) ; level 2 via direct congruence in Z[i]/q^2.
import random
from lib import *
rng = random.Random(3)
ps = split_primes(5, 200)
bad = 0; tot = 0
for trial in range(60):
    P = sorted(rng.sample(ps, 5)); pis = [gauss_prime(p) for p in P]
    for q in [q for q in primerange(3, 120) if q not in P]:
        T = Torus(q); m = T.m
        L = [T.L(pi) for pi in pis]
        for v in monos(5, 3)[::7]:
            z = Zv(pis, v)
            x = sum(a * b for a, b in zip(v, L)) % m
            for s in range(4):
                coll = (x == (s * T.iota) % m)
                div = (Ds(z, s) % q == 0)
                tot += 1; bad += (coll != div)
                # level 2: Z = i^s conj Z mod q^2  <=> q^2 | D_s
                Q2 = q * q
                X, Y = z; is_ = [(1, 0), (0, 1), (-1, 0), (0, -1)][s]
                w = gmul(is_, gconj(z))
                c2 = ((X - w[0]) % Q2 == 0 and (Y - w[1]) % Q2 == 0)
                bad += (c2 != (Ds(z, s) % Q2 == 0)); tot += 1
print("checks", tot, "mismatches", bad)
