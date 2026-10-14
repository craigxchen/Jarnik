# Referee: independent exact check of m_max (T=1), thresholds m_max >= 4w, closed forms of Thm 2.1(b),(e).
import math, sys
from mpmath import mp, mpf, log, sqrt, ln, iv
mp.dps = 50
N_MAX = int(sys.argv[1]) if len(sys.argv) > 1 else 12000

def mmax_exact(d, guess):
    # largest m with 8*C(m+d-1,d) <= 3^d ; monotone in m
    P = 3**d
    ok = lambda m: 8*math.comb(m+d-1, d) <= P
    m = max(1, guess)
    if not ok(m):
        while not ok(m): m -= 1
        return m
    while ok(m+1): m += 1
    return m

def ge4w(m, d, lam):
    # m >= 4 sqrt(8 d lam)  <=> m^2 >= 128 d lam  (m>0); lam = ln(4r) evaluated at high precision
    return mpf(m)**2 >= 128*d*lam

cstar = mpf('0.6417'); kappa = ln(8)/ln(1+1/cstar)
viol_b = 0; viol_ub = 0; last_fail = {'n':None,'nlnn':None,'n2':None}
mm = {}
prev = 1
for n in range(3, N_MAX+1):
    d = n-1
    m = mmax_exact(d, prev if prev > 1 else 1)
    prev = m + 1
    mm[n] = m
    if not (m > cstar*d - kappa): viol_b += 1
    c0 = mpf('0.6417941')  # slightly above root; upper bound check uses c0 root, verify loosely
    if d >= 3 and not (m <= 2*d): viol_ub += 1
    for key, r in (('n', mpf(n)), ('nlnn', n*ln(n)), ('n2', mpf(n)**2)):
        if not ge4w(m, d, ln(4*r)): last_fail[key] = n
print("N_MAX", N_MAX, "violations of m_max > c*d - kappa:", viol_b, " of m_max <= 2d:", viol_ub)
print("last n failing m_max >= 4w:", last_fail)
for n in (1000, 2899, 2900, 3000, 3622, 3623, 4096, 5807, 5808, 5809, 8000, 10000, 12000):
    if n <= N_MAX:
        print(n, "m_max =", mm[n], " m/n = %.6f" % (mm[n]/n))
# kappa, interval constants
iv.dps = 40
g = lambda c: (1+c)*iv.log(1+c) - c*iv.log(c)
c = iv.mpf('0.6417')
print("g(c*) - ln3 in", g(c) - iv.log(3))
print("kappa =", kappa)
t = iv.mpf(316)
print("phi(316) =", c*t - 4*iv.sqrt(8*t), " ln4*phi(316) =", iv.log(4)*(c*t - 4*iv.sqrt(8*t)))
print("6 ln3/sqrt8 =", 6*iv.log(3)/iv.sqrt(8))
# K=316 thresholds: smallest n0 such that n-1 >= 316 ln(4n) for all n>=n0 (increasing beyond 316)
def first(f, start):
    n = start
    while not f(n): n += 1
    return n
print("316 ln 4n threshold:", first(lambda n: n-1 >= 316*ln(4*n), 400), "  316 ln 4n^2 threshold:", first(lambda n: n-1 >= 316*ln(4*mpf(n)**2), 700))
# c0 root
from mpmath import findroot
c0 = findroot(lambda c: (1+c)*log(1+c) - c*log(c) - log(3), 0.64)
print("c0 =", c0)
