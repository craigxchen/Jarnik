"""Prop 4.x: exceptional strata of the density form AC_T (uniform measure on the full ternary cube T).
Rigorous ingredients: sub-Gaussian tail Pr(|S_N| >= t) <= 2 exp(-3 t^2/(4N)) for S_N = sum of N iid uniform{-1,0,1}
(mgf (1+2cosh s)/3 <= e^{s^2/3}, termwise); union bound; Lagrange interpolation degree 2c per box coordinate.
(1) partition strata: q equal classes (independent block sums): mu(Y) >= prod_i (1 - tail_i).
(2) bit strata: q 'bits' with |X_i| = n/2 (any n divisible by 2^q, classes of size n/2^q): mu(Y) >= 1 - q tail.
For (2) also: failure of AC_T at points of {Gamma <= abar} (nearstrata Prop 5.1-type point), W = W_min and W = 3 n ln n.
"""
import math, sys
from common import *

def tail_bound(N, c):
    # Pr(|S_N| > c) = Pr(|S_N| >= c+1)
    return 2 * exp(-3 * mpf(c + 1) ** 2 / (4 * N))

def best_c(N, target):
    # smallest integer c with tail_bound(N,c) <= target
    c = int(sqrt(4 * N * log(2 / mpf(target)) / 3))
    while c > 0 and tail_bound(N, c - 1) <= target:
        c -= 1
    while tail_bound(N, c) > target:
        c += 1
    return c

def partition_qmax(n, r):
    d, w, j, m = wjm(n, r)
    best = 0
    for q in range(1, 400):
        N = d // q  # q classes of size >= floor(d/q) (rest of coordinates go to the base class, value 0)
        # equal targets: each block tail <= 1 - 2^{-1/q}
        tgt = 1 - mpf(2) ** (-mpf(1) / q)
        c = best_c(N, tgt)
        rho = 2 * c * q
        if rho <= m - 1:
            best = q
            bc = (q, N, c, rho)
    return d, m, best, bc

def gamma_point(n, q, eta):
    """Gamma(hat u) upper bound for u = y + eta v, y = binary fraction on 2^q classes of size s = n/2^q."""
    K = 2 ** q
    s = n // K
    N2 = mpf(n) * (n - 1) / 2
    within_pairs = K * mpf(s) * (s - 1) / 2
    # sum_{t=1}^{s-1} (s-t) ln t
    S = mpf(0)
    for t in range(1, s):
        S += (s - t) * log(t)
    within = K * (mpf(s) * (s - 1) / 2 * log(s * (1 + eta) / eta) - S)
    between = (N2 - within_pairs) * log(2 * (1 + eta) * (K - 1))
    return 2 * (within + between) / N2

def bits(n, r, Wlist):
    d, w, j, m = wjm(n, r)
    out = []
    for q in range(1, 40):
        K = 2 ** q
        if n % K or n // K < 2:
            break
        c = best_c(n // 2, mpf(1) / (2 * q))  # mu(Y) >= 1 - q * tail >= 1/2
        rho = 2 * c * q
        exc = rho <= m - 1
        o = m - rho
        res = [q, c, rho, exc, o, K * o > d]
        for W in Wlist:
            if not exc:
                res.append('-'); continue
            le = ln_inv_eps0(n, m, W)  # ln(1/eps0)
            leta = -(le + m * log(4)) / o   # eta_* with 4^m eta^o = eps0
            eta = exp(leta)
            if eta > mpf(1) / (2 * (K - 1)):
                res.append('eta>gap'); continue
            g = gamma_point(n, q, eta)
            abar = W / (2 * d)
            res.append('FAIL(G=%s<=abar=%s)' % (mp.nstr(g, 5), mp.nstr(abar, 5)) if g <= abar else 'ok(G=%s>abar=%s)' % (mp.nstr(g, 5), mp.nstr(abar, 5)))
        out.append(res)
    return d, m, out

if __name__ == '__main__':
    print("== (1) partition strata: largest q with equal classes exceptional for AC_T (rigorous sub-Gaussian) ==")
    for n, r in [(2900, 2900), (4096, 4095), (10**4, 10**4), (10**4, 10**8), (65536, 65535), (10**5, 10**5), (10**6, 10**6), (10**8, 10**8)]:
        d, m, q, bc = partition_qmax(n, r)
        print(n, r, "m'=", m, "q_max>=", q, "(q,N,c,rho)=", bc)
    sys.stdout.flush()
    print("== (2) bit strata, n = 2^k, r = n-1; columns: q c rho exceptional o K*o>d | at W_min | at 3n ln n ==")
    for k in range(12, 31):
        n = 2 ** k
        if k > 22:
            # gamma_point loop over s too long; skip the Gamma test, print exceptional data only
            d, w, j, m = wjm(n, n - 1)
            for q in range(1, 30):
                c = best_c(n // 2, mpf(1) / (2 * q)); rho = 2 * c * q
                if rho <= m - 1:
                    print(k, q, c, rho, True, m - rho, 2 ** q * (m - rho) > d)
            sys.stdout.flush()
            continue
        d, m, out = bits(n, n - 1, [Wmin(n), 3 * n * log(n)])
        for row in out:
            if row[3]:
                print(k, *row)
        sys.stdout.flush()
