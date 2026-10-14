"""Prop 5.3 of nearstrata.md: at Sylvester orders n = 2^k the span of the k columns 1[x_i = 0] contains the point
y_x = sum_i 2^{-i} [x_i = 0], whose values are the n numbers t/n (t = 0..n-1). We compute Gamma(hat u) for it:
pairs at value distance (b-a)/n, |u|_inf = max_x |y_x - y_base| <= (n-1)/n, so
Gamma(hat u) = (2/N) sum_{a<b} ln^+(|u|_inf n /(b-a)) <= (2/N) sum_{t=1}^{n-1} (n-t) ln((n-1)/t)   (base point chosen at t = 0).
Also: exact j, m' = ceil(4w), 2jk and m_max (exact big-integer binomial test) at k = 16, 17, r = n - 1;
and the lower bound 1.5 ln n - 7 for abar of a (hypothetical) cluster with this profile and C <= 1.
"""
import math, mpmath
mpmath.mp.dps = 30

def gamma_binary(k):
    n = 2 ** k; N = n * (n - 1) // 2
    s = mpmath.fsum((n - t) * mpmath.log(mpmath.mpf(n - 1) / t) for t in range(1, n))
    return 2 * s / N

def m_max(n):
    d = n - 1
    rhs = d * mpmath.log(3) - mpmath.log(8)
    lo, hi = 1, 2 * n
    while lo < hi:
        mid = (lo + hi + 1) // 2
        lb = mpmath.loggamma(mid + d) - mpmath.loggamma(d + 1) - mpmath.loggamma(mid)
        if lb <= rhs: lo = mid
        else: hi = mid - 1
    m = lo
    ok = lambda m: 8 * math.comb(m + d - 1, d) <= 3 ** d
    while not ok(m): m -= 1
    while ok(m + 1): m += 1
    return m

for k in (10, 12, 14, 16, 17):
    n = 2 ** k
    print(f"k={k}: Gamma(binary point) <= {float(gamma_binary(k)):.4f};  1.5 ln n - 7 = {1.5*math.log(n)-7:.3f}")
for k in (15, 16, 17):
    n = 2 ** k; d = n - 1; r = n - 1
    w = math.sqrt(8 * d * math.log(4 * r)); j = math.floor(w / 2); mp = math.ceil(4 * w)
    mm = m_max(n)
    print(f"k={k}: n={n}, w={w:.3f}, j={j}, ceil(4w)={mp}, 2jk={2*j*k}, m_max={mm}, window (2jk, m_max] nonempty: {2*j*k < mm}, ceil(4w) <= 2jk: {mp <= 2*j*k}")
