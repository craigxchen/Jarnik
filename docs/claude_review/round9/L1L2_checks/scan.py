"""Enumerate normalised actual clusters (M >= Mmin points on an arc of length <= Cmax sqrt R).

Usage: python3 scan.py Nmax T Cmax Mmin out.pkl
  N runs over split-only integers (all prime factors = 1 mod 4) <= Nmax with tau(N) >= T.
  A cluster is a set of angularly consecutive lattice points of x^2+y^2=N with
  (theta_max - theta_min) N^{1/4} <= Cmax; only clusters whose points have gcd 1 (normalised)
  are kept; clusters are deduplicated under the 8 symmetries of Z[i].
Exact data stored per cluster: N, factorisation, the points.  Statistics are computed in stats.py.
"""
import sys, pickle
from math import atan2, pi
from glib import primes_upto, points_of_norm, ggcd, norm

Nmax = int(float(sys.argv[1])); T = int(sys.argv[2]); Cmax = float(sys.argv[3]); Mmin = int(sys.argv[4])
out = sys.argv[5]

P1 = [p for p in primes_upto(max(10 ** 6, Nmax // (5 * 13 * 17 if T >= 16 else 65) + 1)) if p % 4 == 1]


def gen(idx, cur, fac, tau):
    # yield (N, fac) with tau(N) >= T; prune with tau * 2^floor(log(budget)/log p) < T
    from math import log as _log
    if tau >= T:
        yield cur, dict(fac)
    for k in range(idx, len(P1)):
        p = P1[k]
        if cur * p > Nmax:
            break
        budget = Nmax // cur
        if tau * 2 ** int(_log(budget) / _log(p) + 1e-9) < T:
            break
        e = 0; q = cur
        while q * p <= Nmax:
            q *= p; e += 1
            fac[p] = e
            yield from gen(k + 1, q, fac, tau * (e + 1))
            del fac[p]


def canon(pts):
    best = None
    for c in (False, True):
        for r in range(4):
            Q = []
            for z in pts:
                w = (z[0], -z[1]) if c else z
                for _ in range(r):
                    w = (-w[1], w[0])
                Q.append(w)
            key = tuple(sorted(Q))
            if best is None or key < best:
                best = key
    return best


seen = {}
nN = 0
for N, fac in gen(0, 1, {}, 1):
    nN += 1
    pts = points_of_norm(fac)
    n = len(pts)
    ang = [atan2(z[1], z[0]) % (2 * pi) for z in pts]
    s = N ** 0.25
    for i in range(n):
        # extend window
        j = i
        while True:
            jn = (j + 1) % n
            if jn == i:
                break
            d = (ang[jn] - ang[i]) % (2 * pi)
            if d * s > Cmax:
                break
            j = jn
        # all sub-windows [i, j'] with >= Mmin points starting at i
        cnt = (j - i) % n + 1
        for m in range(Mmin, cnt + 1):
            cl = [pts[(i + t) % n] for t in range(m)]
            g = cl[0]
            for z in cl[1:]:
                g = ggcd(g, z)
            if norm(g) != 1:
                continue
            Cc = ((ang[(i + m - 1) % n] - ang[i]) % (2 * pi)) * s
            key = canon(cl)
            if key not in seen:
                seen[key] = (N, fac, cl, Cc)
print("N scanned", nN, "clusters", len(seen))
with open(out, "wb") as f:
    pickle.dump(list(seen.values()), f)
