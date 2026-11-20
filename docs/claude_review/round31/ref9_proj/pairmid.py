"""Check T(e_i+e_j)*{1,2} = |prim(z_i+z_j)|^2 on all pairs of normalised clusters (z_i != -z_j)."""
import sys; sys.path.insert(0, '/tmp/claude-0/-home-user-Jarnik/ba4929a4-8665-5790-985a-3f09130a14b9/scratchpad/round31/ref9_proj')
from structured import *
from math import gcd
cnt = 0; bad = 0
for N in list(range(5, 30000)) + [1176852625, 235370525, 537606725]:
    f = factorint(N)
    if any(p % 4 == 3 and e % 2 for p, e in f.items()): continue
    pts = all_points(N)
    if len(pts) < 3 or len(pts) > 200: continue
    # treat consecutive triples as clusters
    for t in range(len(pts)):
        cl = [pts[t], pts[(t+1) % len(pts)], pts[(t+2) % len(pts)]]
        g, zs, Np = normalise(cl)
        if Np == 1: continue
        primes, pis, e, rows = rows_of(zs, Np)
        for i, j in combinations(range(3), 2):
            s = (zs[i][0]+zs[j][0], zs[i][1]+zs[j][1])
            if s == (0, 0): continue
            G = gcd(*s); h = nrm(s)//(G*G)
            c = [0]*3; c[i] += 1; c[j] += 1
            T, b, u, v = certificate(zs, Np, primes, pis, e, rows, c)
            cnt += 1
            if nrm(b) != h: bad += 1
print("pairs", cnt, "with |b|^2 != |primitive(z_i+z_j)|^2:", bad)
