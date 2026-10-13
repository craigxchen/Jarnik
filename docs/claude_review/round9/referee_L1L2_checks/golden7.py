"""Referee: golden family member n=7 (N'=55728887765525): margins via own code (points by Gaussian factorisation)."""
import itertools
from math import log, atan2, pi
from refgauss import *
from records import profile
N = 55728887765525
fac = factor(N); print(fac)
pts = [(1, 0)]
for p, e in fac.items():
    pp = gauss_prime(p); new = []
    for z in pts:
        for a in range(e + 1):
            new.append(gmul(z, gmul(gpow(pp, a), gpow(gconj(pp), e - a))))
    pts = new
pts = [gmul(z, u) for z in pts for u in UNITS]
assert all(gnorm(z) == N for z in pts) and len(set(pts)) == len(pts)
pts.sort(key=lambda z: atan2(z[1], z[0]) % (2 * pi))
ang = [atan2(z[1], z[0]) % (2 * pi) for z in pts]; n = len(pts); s = N ** 0.25
best = min(((ang[(i + 7) % n] - ang[i]) % (2 * pi) * s, i) for i in range(n))
C, i = best
cl = [pts[(i + t) % n] for t in range(8)]
u, G = normalise(cl); Np, f2, pis, rows, ks = profile(u); W = log(Np)
mus = sorted(log(Np) - log(gnorm(ggcd(u[x], u[y]))) - W / 2 for x, y in itertools.combinations(range(8), 2))
mbar = sum(mus) / len(mus)
print(f"C={C:.4f} |G|^2={gnorm(G)} W={W:.3f} W/2={W/2:.3f} mubar={mbar:.3f} max={mus[-1]:.3f} K={mus[-1]/mbar:.3f}")
print("sorted margins:", [round(m, 2) for m in mus])
