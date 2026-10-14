"""Round 12, T = 1 parameters (exact where stated).

m_max(n): largest m with 8*binom(m+n-2, n-1) <= 3^(n-1)  (exact big-integer check at the boundary).
Closed form (Lemma 2.2 of t1rigour.md): m_max >= c* d - kappa, d = n-1, c* = 0.6417, kappa = ln 8 / ln(1 + 1/c*).
w = sqrt(8 d ln 4r).  Thresholds for m_max >= 4w; closed-form sufficient condition d >= 318 ln 4r.
W_1 exact form and closed form for m' = m_max and m' = ceil(4w).
"""
import math, sys
import mpmath
mpmath.mp.dps = 50

def g(c):
    c = mpmath.mpf(c)
    return (1 + c) * mpmath.log(1 + c) - c * mpmath.log(c)

# c0 and c*
c0 = mpmath.findroot(lambda c: g(c) - mpmath.log(3), 0.64)
cstar = mpmath.mpf('0.6417')
kappa = mpmath.log(8) / mpmath.log(1 + 1 / cstar)
print("c0 =", mpmath.nstr(c0, 20), " g(c*) - ln3 =", mpmath.nstr(g(cstar) - mpmath.log(3), 10),
      " kappa = ln8/ln(1+1/c*) =", mpmath.nstr(kappa, 10))
assert g(cstar) < mpmath.log(3)

def ok(m, n):
    return 8 * math.comb(m + n - 2, n - 1) <= 3 ** (n - 1)

def m_max(n):
    d = n - 1
    rhs = d * mpmath.log(3) - mpmath.log(8)
    lo, hi = 1, 2 * n
    while lo < hi:
        mid = (lo + hi + 1) // 2
        lb = mpmath.loggamma(mid + n - 1) - mpmath.loggamma(n) - mpmath.loggamma(mid)
        if lb <= rhs:
            lo = mid
        else:
            hi = mid - 1
    m = lo
    # exact confirmation
    while not ok(m, n):
        m -= 1
    while ok(m + 1, n):
        m += 1
    return m

N0, N1 = int(sys.argv[1]) if len(sys.argv) > 1 else 1000, int(sys.argv[2]) if len(sys.argv) > 2 else 9000
rs = {'r=n': lambda n: n, 'r=n ln n': lambda n: n * math.log(n), 'r=n^2': lambda n: n * n}
fail_last = {k: None for k in rs}
cf_fail = 0
minratio = 10.0
rows = []
for n in range(N0, N1 + 1):
    d = n - 1
    m = m_max(n)
    cf = float(cstar) * d - float(kappa)
    if m < cf:
        cf_fail += 1
    minratio = min(minratio, m / n)
    for k, f in rs.items():
        w = math.sqrt(8 * d * math.log(4 * f(n)))
        if m < 4 * w:
            fail_last[k] = n
    if n % 500 == 0 or n in (2900, 3000, 3623, 5807, 5808, 5809):
        rows.append((n, m, m / n, m - cf))
print(f"n in [{N0},{N1}]: closed-form bound m_max >= c* d - kappa violated {cf_fail} times; min m/n = {minratio:.5f}")
for k in rs:
    print(f"  {k}: last n with m_max < 4w in range: {fail_last[k]}")
print("n, m_max, m_max/n, m_max - (c* d - kappa)")
for r in rows:
    print(f"  {r[0]}, {r[1]}, {r[2]:.5f}, {r[3]:.3f}")
