"""Detail for the 9-point record run: exceptional pairs of Thm 3.2 and their margins; count of N in tangent.py."""
import itertools, random
from math import log
from glib import ggcd, norm, sub, primes_upto
import tangent
N = 1176852625
runs = tangent.all_runs(N, 8)
cl = max(runs, key=len)
cl, Nn = tangent.normalise(cl)
W = log(Nn)
out = []
tangent.analyse(cl, "rec", out)
r = out[0]
v = r[6]; iv = (-v[1], v[0])
print("N'", Nn, "M", len(cl), "h*", r[5], "v*", v)
for x, y in itertools.combinations(range(len(cl)), 2):
    u = tangent.prim(sub(cl[x], cl[y]))
    if u[0] * iv[1] - u[1] * iv[0] == 0:
        g = ggcd(cl[x], cl[y])
        print("exceptional", x, y, "mu - W/2 =", round(W / 2 - log(norm(g)) - W / 2, 6), "Norm g", norm(g))
# reproduce the N count of tangent.py main
split = [p for p in primes_upto(400) if p % 4 == 1]
Ns = set(); random.seed(1)
for _ in range(220):
    k = random.randint(3, 6); ps = random.sample(split[:30], k); n = 1
    for p in ps: n *= p ** random.randint(1, 2)
    if n < 10 ** 16: Ns.add(n)
print("random N count", len(Ns))
