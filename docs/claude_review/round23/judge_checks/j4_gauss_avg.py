# Judge J4 (unrefereed): exact check of the Gaussian random-5-subset averaging on the record's 10-point cluster
# n = R^2 = 1176852625 = 5^3 13^2 17 29 113. For each 5-subset S: Gaussian gcd g_S, reduced N''=n/|g_S|^2,
# identity N'' = prod_p p^(max_S a_p - min_S a_p), and mean of log N''/log n vs 15/16 + 15/(16(M-4)).
from math import isqrt, atan2, log, sqrt
from itertools import combinations
from fractions import Fraction as Fr
n = 1176852625
pts = []
for x in range(0, isqrt(n)+1):
    y2 = n - x*x; y = isqrt(y2)
    if y*y == y2:
        for sx in (1,-1):
            for sy in (1,-1):
                pts.append((sx*x, sy*y))
pts = sorted(set(pts), key=lambda p: atan2(p[1], p[0]))
R = sqrt(n); k = 10
best = None
L = len(pts)
for i in range(L):
    cl = [pts[(i+j) % L] for j in range(k)]
    span = sqrt((cl[0][0]-cl[-1][0])**2 + (cl[0][1]-cl[-1][1])**2)
    if best is None or span < best[0]: best = (span, cl)
span, cl = best
print("lattice points:", L, " best 10-point span / sqrt(R) =", span/sqrt(R))
def gmod(a, b):  # remainder of Gaussian integers a mod b
    (ar, ai), (br, bi) = a, b
    nb = br*br + bi*bi
    qr = ar*br + ai*bi; qi = ai*br - ar*bi          # a * conj(b)
    q = (round(Fr(qr, nb)), round(Fr(qi, nb)))
    return (ar - (q[0]*br - q[1]*bi), ai - (q[0]*bi + q[1]*br))
def ggcd(a, b):
    while b != (0, 0): a, b = b, gmod(a, b)
    return a
def norm(a): return a[0]**2 + a[1]**2
# split primes and a Gaussian prime pi for each
split = {5: (2, 1), 13: (3, 2), 17: (4, 1), 29: (5, 2), 113: (8, 7)}
def val(z, pi):
    c = 0
    while True:
        pr, pi_ = pi; nb = norm(pi)
        qr = z[0]*pr + z[1]*pi_; qi = z[1]*pr - z[0]*pi_
        if qr % nb or qi % nb: return c
        z = (qr // nb, qi // nb); c += 1
M = len(cl); ratios = []; ok = True
for S in combinations(range(M), 5):
    g = (0, 0)
    for j in S: g = ggcd(g, cl[j])
    N2 = n // norm(g); assert n % norm(g) == 0
    prod = 1
    for p, pi in split.items():
        a = [val(cl[j], pi) for j in S]
        prod *= p ** (max(a) - min(a))
    ok &= (prod == N2)
    ratios.append(log(N2) / log(n))
mean = sum(ratios) / len(ratios)
print("identity N''=prod p^(max a-min a) on all", len(ratios), "subsets:", ok)
print("mean log N''/log n = %.6f ; bound 15/16+15/(16(M-4)) = %.6f ; max = %.6f ; min = %.6f" %
      (mean, 15/16 + 15/(16*(M-4)), max(ratios), min(ratios)))
