# Referee: independent exact check of the good-prime dictionary (D5) on actual primitive clusters.
# For p not dividing 2L, alpha_j=v_pi(X_j+iL), gamma_j=v_pibar(X_j+iL) (j finite; anchor alpha=gamma=0):
#   min(alpha_j,gamma_j)=0;  v_p(N)=max alpha+max gamma;  v_p(X_j^2+L^2)=alpha_j+gamma_j;
#   v_p(X_i-X_j)=min(alpha_i,alpha_j)+min(gamma_i,gamma_j);
#   a_jp - a_0p = alpha_j - gamma_j (up to the global choice pi <-> pibar), a_0p = max gamma.
import sys, math, random
from fractions import Fraction
from math import gcd
sys.path.insert(0, '.')
from rlib import *

def vp_int(n, p):
    if n == 0: return 10 ** 9
    v = 0
    while n % p == 0: n //= p; v += 1
    return v
random.seed(int(sys.argv[1]) if len(sys.argv) > 1 else 4)
SPLIT = [5, 13, 17, 29, 37, 41, 53, 61]
tested = 0; checks = 0; fails = 0
for trial in range(2500):
    ps = random.sample(SPLIT, random.randint(1, 4)); N0 = 1
    for p in ps: N0 *= p ** random.randint(1, 3)
    if N0 > 10 ** 8: continue
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
    X = [cot[0, j] * L for j in range(1, M)]; X = [int(x) for x in X]
    tested += 1
    for p, e in factor(N).items():
        if (2 * L) % p == 0: continue
        pi = gauss_prime_above(p); pib = gconj(pi)
        al = [vpi((x, L), pi) for x in X]; ga = [vpi((x, L), pib) for x in X]
        a = [vpi(z, pi) for z in P]
        conds = [all(min(u, w) == 0 for u, w in zip(al, ga)),
                 e == max(al + [0]) + max(ga + [0]),
                 all(vp_int(x * x + L * L, p) == u + w for x, u, w in zip(X, al, ga)),
                 all(vp_int(X[i] - X[j], p) == min(al[i], al[j]) + min(ga[i], ga[j])
                     for i in range(len(X)) for j in range(i + 1, len(X)))]
        o1 = all(a[j + 1] - a[0] == al[j] - ga[j] for j in range(len(X))) and a[0] == max(ga + [0])
        o2 = all(a[j + 1] - a[0] == ga[j] - al[j] for j in range(len(X))) and a[0] == max(al + [0])
        conds.append(o1 or o2)
        checks += 1
        if not all(conds): fails += 1; print('D5 FAIL', N, p, X, L, conds)
print('clusters', tested, 'good-prime checks', checks, 'failures', fails)
