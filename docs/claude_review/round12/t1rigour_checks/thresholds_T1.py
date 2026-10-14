"""Closed-form thresholds for m_max >= 4w at T = 1 (Lemma 2.2/2.3 of t1rigour.md) and Theorem 3.1 constants W_1.
Closed form used: m_max >= c* d - kappa (c* = 0.6417, kappa = 2.2137), w = sqrt(8 d ln 4r).
f(d, lam) = c* d - kappa - 4 sqrt(8 d lam) >= 0 is sufficient for m_max >= 4w.
"""
import math
import mpmath
mpmath.mp.dps = 40
cs, kap = 0.6417, math.log(8) / math.log(1 + 1 / 0.6417)

def f(d, lam):
    return cs * d - kap - 4 * math.sqrt(8 * d * lam)

# K-criterion
for K in (315, 316, 318, 320):
    print(f"K={K}: c*K - 4 sqrt(8K) = {cs*K - 4*math.sqrt(8*K):.4f}; needs lam >= {kap/(cs*K - 4*math.sqrt(8*K)) if cs*K > 4*math.sqrt(8*K) else float('inf'):.4f} (ln 4 = {math.log(4):.4f})")

def first_n(lamf, start=10):
    # smallest n such that f(d, lam(n)) >= 0 for this and all larger n (monotone beyond; checked to 10^7 on a grid)
    n = start
    while f(n - 1, lamf(n)) < 0:
        n += 1
    # check monotone holding on [n, 2*10^5] exactly and on a log grid to 1e9
    for k in range(n, 200001):
        assert f(k - 1, lamf(k)) >= 0, k
    x = 2e5
    while x < 1e9:
        assert f(x - 1, lamf(x)) >= 0
        x *= 1.01
    return n

print("closed-form n_0 (sufficient):")
print("  r = n     :", first_n(lambda n: math.log(4 * n)))
print("  r = n ln n:", first_n(lambda n: math.log(4 * n * math.log(n))))
print("  r = n^2   :", first_n(lambda n: math.log(4 * n * n)))
print("  r = n via d >= 316 ln 4r:", next(n for n in range(10, 10**5) if n - 1 >= 316 * math.log(4 * n)))
print("  r = n^2 via d >= 316 ln 4r:", next(n for n in range(10, 10**5) if n - 1 >= 316 * math.log(4 * n * n)))

def m_max(n):
    d = n - 1
    rhs = d * mpmath.log(3) - mpmath.log(8)
    lo, hi = 1, 2 * n
    while lo < hi:
        mid = (lo + hi + 1) // 2
        if mpmath.loggamma(mid + n - 1) - mpmath.loggamma(n) - mpmath.loggamma(mid) <= rhs:
            lo = mid
        else:
            hi = mid - 1
    return lo

print("\nTheorem 3.1 constants (r = n, C = 1): W1ex(m') = [6 d ln3 + 2 m' ln d - 2 ln m'!]/(m'/2 - w);"
      " W1cf(m') = 24 d ln3/m' + 8 ln(e d/m')")
print("n, w, m'=ceil(4w), m_max, W1ex(ceil4w), W1cf(ceil4w), W1ex(m_max), W1cf(m_max), 3n ln n - 14n")
for n in (2900, 3000, 5000, 10**4, 10**5, 10**6, 10**7):
    d = n - 1
    w = math.sqrt(8 * d * math.log(4 * n))
    mp_ = math.ceil(4 * w)
    mm = m_max(n)
    def W1ex(m):
        return (6 * d * math.log(3) + 2 * m * math.log(d) - 2 * math.lgamma(m + 1)) / (m / 2 - w)
    def W1cf(m):
        return 24 * d * math.log(3) / m + 8 * math.log(math.e * d / m)
    print(f"{n}, {w:.1f}, {mp_}, {mm}, {W1ex(mp_):.2f}, {W1cf(mp_):.2f}, {W1ex(mm):.2f}, {W1cf(mm):.2f}, {3*n*math.log(n)-14*n:.4g}")
