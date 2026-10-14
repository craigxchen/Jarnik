"""Recheck of J2 numbers (round11/README.md §5.2): W_TC = 2(m-1) ln 2 + 2 ln((D/2)/(m+1)) at T = n and T = 1."""
import math, mpmath
mpmath.mp.dps = 50
def m_of(n, T):
    d = n - 1
    rhs = d * mpmath.log(2 * T + 1) - mpmath.log(8)
    lo, hi = 1, 4 * n * T + 10
    while lo < hi:
        mid = (lo + hi + 1) // 2
        if mpmath.loggamma(mid + n - 1) - mpmath.loggamma(n) - mpmath.loggamma(mid) <= rhs:
            lo = mid
        else:
            hi = mid - 1
    m = lo
    if n <= 4000:
        R = (2 * T + 1) ** (n - 1)
        assert 8 * math.comb(m + n - 2, n - 1) <= R < 8 * math.comb(m + n - 1, n - 1)
    return m
for n in (2137, 3000):
    m = m_of(n, n)
    Dh = (n - 1) * n
    WTC = 2 * (m - 1) * math.log(2) + 2 * math.log(Dh / (m + 1))
    Wmin = 3 * n * math.log(n) - 14 * n
    print(f"T=n, n={n}: m={m}, m/(n(n-1))={m/(n*(n-1)):.4f}, W_TC={WTC:.4g}, W_TC/n^2={WTC/n**2:.4f}, W_min={Wmin:.5g}, ratio={WTC/Wmin:.1f}")
for n in (3000, 4000, 8000):
    m = m_of(n, 1)
    WTC = 2 * (m - 1) * math.log(2) + 2 * math.log((n - 1) / (m + 1))
    Wmin = 3 * n * math.log(n) - 14 * n
    room = Wmin / 4 + (1 - m) * math.log(2) - math.log((n - 1) / (m + 1))
    print(f"T=1, n={n}: m={m}, W_TC={WTC:.2f} (=1.39m? {WTC/m:.4f}), W_min={Wmin:.1f}, log room at W_min={room:.1f}")
