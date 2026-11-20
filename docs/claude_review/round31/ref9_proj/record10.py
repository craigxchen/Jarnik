import sys; sys.path.insert(0, '/tmp/claude-0/-home-user-Jarnik/ba4929a4-8665-5790-985a-3f09130a14b9/scratchpad/round31/ref9_proj')
from cert import *
from math import isqrt
N = 1176852625
print(factorint(N))
# enumerate all lattice points of N and find the 10-point cluster near the diagonal
pts = []
for x in range(isqrt(N)+1):
    y2 = N - x*x; y = isqrt(y2)
    if y*y == y2:
        for sx in (1, -1):
            for sy in (1, -1):
                pts.append((sx*x, sy*y))
pts = sorted(set(pts), key=lambda z: atan2(z[1], z[0]))
R = sqrt(N)
# points within angle pi/4 +- 0.01
cl = [z for z in pts if abs(atan2(z[1], z[0]) - pi/4) < 0.03]
print(len(cl), cl)
a, th, wid = angles(cl)
print("arc length / sqrt R =", R*wid/sqrt(R), " chord C' =", max(sqrt(nrm((z[0]-w[0], z[1]-w[1]))) for z in cl for w in cl)/N**0.25)
g, zs, Np = normalise(cl)
print("gcd", g, "N'", Np)
primes, pis, e, rows = rows_of(zs, Np)
print(primes, e); [print(r) for r in rows]
M = len(zs); T1 = []
for i, j in combinations(range(M), 2):
    c = [0]*M; c[i] += 1; c[j] += 1
    T, b, u, v = certificate(zs, Np, primes, pis, e, rows, c)
    xs = projections(zs, b)
    if T == 1:
        T1.append((i, j, b, u, xs))
print("pairs with T=1:", len(T1))
for t in T1: print(t)
# all integers in the interval with 2N|b|^2 - x^2 square
b = T1[0][2]; xs = T1[0][4]; D = Np*nrm(b)
cand = [x for x in range(xs[0], xs[-1]+1) if (D-x*x) >= 0 and isqrt(D-x*x)**2 == D-x*x]
print("all norm lifts in [min,max]:", cand)
print("sieve terms (lhs, sum log|diff|, n(n-1)logX):", sieve_terms(xs, D, M))
print("(2+X)/M^1.5 =", (2+xs[-1]-xs[0])/M**1.5, " M^2 exp(-bad) =", M**2*2.718281828**(-bad_sum(Np, M)))
