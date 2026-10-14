"""Counting lower bound for rho(Y_t(j)) as j -> infinity (refJ §6.4 question).
rho(P) >= min{e : #{c in [0,2j]^b : |c| <= e} >= |P|}  (S/J(P) is a quotient of S/(y^{2j+1}), h_P(rho) = |P|).
|Y_t(j)| = |Lambda_t cap [-j,j]^b| >= (2^{(t-1)b}/I_t) * floor((2j+1)/2^{t-1})^b  (Lambda_t contains 2^{t-1} Z^b),
I_t = [Z^b : Lambda_t] = 2^{t-1+(t-2)(2^{t-1}-1)}  (det of the PG(t-1,2) incidence matrix; SNF check t<=4 in
sylvester_lattices_out.txt).  As j -> oo: rho/j >= 2 kappa_t with IrwinHall_b(kappa_t) = 1/I_t (exact rational CDF).
Also prints finite-j bounds for j = 50, 100, 1000 using exact integer counts."""
from fractions import Fraction
from math import comb, floor
import mpmath as mp
mp.mp.dps = 50

def I_t(t):
    return 2 ** (t - 1 + (t - 2) * (2 ** (t - 1) - 1))

def irwin_hall_cdf(b, x):
    # P(U_1+..+U_b <= x)
    s = mp.mpf(0)
    for k in range(0, int(mp.floor(x)) + 1):
        s += (-1) ** k * comb(b, k) * (mp.mpf(x) - k) ** b
    return s / mp.factorial(b)

def kappa(t):
    b = 2 ** t - 1
    target = mp.mpf(1) / I_t(t)
    lo, hi = mp.mpf(0), mp.mpf(b) / 2
    for _ in range(200):
        mid = (lo + hi) / 2
        if irwin_hall_cdf(b, mid) >= target:
            hi = mid
        else:
            lo = mid
    return hi

def capped_count(b, cap, e):
    # number of c in [0,cap]^b with sum <= e  (inclusion-exclusion)
    tot = 0
    for k in range(0, b + 1):
        r = e - k * (cap + 1)
        if r < 0:
            break
        tot += (-1) ** k * comb(b, k) * comb(r + b, b)
    return tot

def finite_bound(t, j):
    b = 2 ** t - 1
    M = 2 ** (t - 1)
    Ylow = (M ** b // I_t(t)) * ((2 * j + 1) // M) ** b
    e = 0
    lo, hi = 0, 2 * j * b
    while lo < hi:
        mid = (lo + hi) // 2
        if capped_count(b, 2 * j, mid) >= Ylow:
            hi = mid
        else:
            lo = mid + 1
    return lo

for t in range(2, 8):
    k = kappa(t)
    print(f't={t} b={2**t-1} log2 I_t={t-1+(t-2)*(2**(t-1)-1)} kappa_t={mp.nstr(k,8)}  asymptotic rho_t/j >= {mp.nstr(2*k,8)}',
          '  finite j=50,100,1000:', [finite_bound(t, j) / j for j in (50, 100, 1000)] if t <= 6 else '')
