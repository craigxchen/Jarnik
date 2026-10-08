# Referee: exact check of the chart master identity (M3)
#   N^{Q*} * L^{2P} * prod_{1<=i<j<=k}(X_i-X_j)^2 * prod_e eps_e = Delta_* * prod_e s_e^2 * prod_i (X_i^2+L^2)^k
# (k = M-1 finite coordinates, anchor z_0 at an endpoint, L = L_0), and of the inequality (M4)
#   N^{Q*} >= 2^{-Q*} (A/L)^{2P} Delta_* prod s^2   (exact rational form), on actual primitive clusters.
import sys, math, random
from fractions import Fraction
from math import gcd
sys.path.insert(0, '.')
from rlib import *

random.seed(int(sys.argv[1]) if len(sys.argv) > 1 else 2)
SPLIT = [5, 13, 17, 29, 37, 41, 53, 61, 73]
ok = 0; ok4 = 0; tested = 0
for trial in range(3000):
    ps = random.sample(SPLIT, random.randint(1, 4)); N0 = 1
    for p in ps: N0 *= p ** random.randint(1, 3)
    if N0 > 2 * 10 ** 8: continue
    pts = sorted(points_on_circle(N0), key=angle); n = len(pts)
    M = random.randint(3, min(8, n // 2))
    st = random.randrange(n); W = [pts[(st + t) % n] for t in range(M)]
    if (angle(W[-1]) - angle(W[0])) % (2 * math.pi) >= math.pi - 1e-9: continue
    g = ggcd_list(W); P = [gdivmod_exact(z, g) for z in W]; N = gnorm(P[0])
    if N == 1: continue
    cot = {}
    for i in range(M):
        for j in range(i + 1, M):
            a, b = P[i], P[j]
            cot[i, j] = Fraction(N + a[0] * b[0] + a[1] * b[1], a[0] * b[1] - a[1] * b[0])
    L = 1
    for c in cot.values(): L = L * c.denominator // gcd(L, c.denominator)
    X = [abs(cot[0, j]) * L for j in range(1, M)]
    assert all(x.denominator == 1 for x in X)
    X = [int(x) for x in X]; k = M - 1; Pn = M * (M - 1) // 2; Qs = M * M // 4
    prod_s2 = 1; prod_eps = 1
    for (i, j), c in cot.items():
        prod_s2 *= c.denominator ** 2
        prod_eps *= 2 if (P[i][0] - P[j][0]) % 2 else 1
    # Delta_* = prod p^{(4 D_p - e_p (M^2 - 4 Q*))/4}
    Dstar = 1
    for p, e in factor(N).items():
        pi = gauss_prime_above(p); lv = [vpi(z, pi) for z in P]
        D4 = sum((2 * sum(1 for v in lv if v < tau) - M) ** 2 for tau in range(1, e + 1))
        ex4 = D4 - e * (M * M - 4 * Qs)
        assert ex4 % 4 == 0 and ex4 >= 0
        Dstar *= p ** (ex4 // 4)
    lhs = N ** Qs * L ** (2 * Pn) * prod_eps
    for i in range(k):
        for j in range(i + 1, k): lhs *= (X[i] - X[j]) ** 2
    rhs = Dstar * prod_s2
    for x in X: rhs *= (x * x + L * L) ** k
    tested += 1
    if lhs == rhs: ok += 1
    else: print('M3 FAIL', N, M, W)
    A = min(X)
    if Fraction(N ** Qs) >= Fraction(A, L) ** (2 * Pn) * Dstar * prod_s2 / 2 ** Qs: ok4 += 1
    else: print('M4 FAIL', N, M, W)
print('clusters', tested, '(M3) exact identity holds on', ok, '; (M4) holds on', ok4)
