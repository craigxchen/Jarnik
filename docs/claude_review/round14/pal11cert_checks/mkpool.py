"""Pool of generators for the product-certificate search on pal11 (L=2, m=8).
usage: python3 mkpool.py maxsupp maxentry out.txt
A generator is a balanced primitive z in Z^12 (row 0 = base), |z_x| <= maxentry for all 12 rows, support 2..maxsupp,
sign-normalised (first nonzero entry > 0); its linear form is <z',u>, z' = z[1:].  Every primitive z' in Z^11 (up to
sign) arises from exactly one such z.  Writes: pool size, test-direction count, then per generator: orbit-rep flag,
z (12 ints), number of touched test directions, their indices.  Orbit = Z_5 (multipliers x -> a x, a square mod 11)."""
import sys, itertools
from math import gcd
import numpy as np
maxsupp, maxentry, out = int(sys.argv[1]), int(sys.argv[2]), sys.argv[3]
dat = np.load('dirs_pal11.npz'); D = dat['D']; tau = dat['tau']
n = 12
pool = []
for s in range(2, maxsupp + 1):
    for S in itertools.combinations(range(n), s):
        for vals in itertools.product([e for e in range(-maxentry, maxentry + 1) if e], repeat=s):
            if sum(vals) != 0 or vals[0] < 0: continue
            g = 0
            for v in vals: g = gcd(g, abs(v))
            if g != 1: continue
            z = [0]*n
            for x, v in zip(S, vals): z[x] = v
            pool.append(tuple(z))
pool.sort(key=lambda z: (sum(1 for x in z if x), sum(abs(x) for x in z), [-abs(x) for x in z], z))
index = {z: i for i, z in enumerate(pool)}
def norm(z):
    f = next(x for x in z if x)
    return tuple(z) if f > 0 else tuple(-x for x in z)
def act(a, z):  # row 1+i -> row 1+(a i mod 11), base fixed
    w = [0]*n; w[0] = z[0]
    for i in range(11): w[1 + (a*i) % 11] = z[1+i]
    return norm(w)
rep = []
for i, z in enumerate(pool):
    orb = [index[act(a, z)] for a in (1, 3, 4, 5, 9)]
    rep.append(1 if i == min(orb) else 0)
P = np.array(pool, dtype=np.int64)[:, 1:]
NZ = (P @ D.T) != 0
with open(out, 'w') as f:
    f.write(f'{len(pool)} {len(D)}\n')
    f.write(' '.join(str(int(t)) for t in tau) + '\n')
    for i, z in enumerate(pool):
        idx = np.nonzero(NZ[i])[0]
        f.write(f'{rep[i]} ' + ' '.join(map(str, z)) + f' {len(idx)} ' + ' '.join(map(str, idx)) + '\n')
print('pool', len(pool), 'orbit reps', sum(rep), 'dirs', len(D))
