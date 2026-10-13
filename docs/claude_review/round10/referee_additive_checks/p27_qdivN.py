"""Referee remark: Prop 2.7(2) also for odd primes q | N (normalised clusters): q | I iff exactly 2 residues mod q."""
import itertools, random, math
import mpmath as mp
from sympy import factorint
from glib_ref import *
from actual_checks import cross
rng = random.Random(5)
recs = [1176852625, 537606725, 235370525]
Ns = [(N, factorint(N)) for N in recs] + [random_split_N(rng, 10 ** 13) for _ in range(150)]
tot = fails = qdiv = 0
for (N, fac) in Ns:
    pts = points_on_circle(fac)
    P = sorted(((z, angle(z)) for z in pts), key=lambda t: t[1])
    for C in (2, 4, 6, 8):
        D = mp.mpf(C) * mp.power(N, mp.mpf(-1) / 4)
        for win in windows(P, D):
            if len(win) < 3:
                continue
            zs, g = normalise([P[k][0] for k in win])
            Np = gnorm(zs[0])
            I = 0
            for (x, y, w) in itertools.combinations(range(len(zs)), 3):
                I = math.gcd(I, abs(cross(gsub(zs[y], zs[x]), gsub(zs[w], zs[x]))))
            for q in factorint(Np):
                tot += 1
                nres = len(set((z[0] % q, z[1] % q) for z in zs))
                if (I % q == 0) != (nres == 2):
                    fails += 1
                qdiv += (I % q == 0)
print('pairs (cluster, q | N) tested', tot, 'with q | I:', qdiv, 'failures of [q | I iff exactly 2 residues]:', fails)
