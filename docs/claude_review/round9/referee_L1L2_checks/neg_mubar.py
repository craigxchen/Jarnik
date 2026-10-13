"""Referee: clusters of the scan with mubar <= 0 (L1 fails for every K); their A = W/(M log M) and C."""
import pickle, itertools
from math import log
from refgauss import *
data = pickle.load(open('ref_1e9.pkl', 'rb'))
out = []
for (N, fac, cl, C) in data:
    M = len(cl); u, G = normalise(cl); Np = gnorm(u[0]); W = log(Np)
    mus = [log(Np) - log(gnorm(ggcd(u[x], u[y]))) - W / 2 for x, y in itertools.combinations(range(M), 2)]
    mb = sum(mus) / len(mus)
    if mb <= 1e-12:
        out.append((W / (M * log(M)), N, M, round(C, 3), round(mb, 4), round(max(mus), 3)))
out.sort()
print(len(out), "clusters with mubar <= 0; smallest A first:")
for r in out[:19]:
    print(r)
