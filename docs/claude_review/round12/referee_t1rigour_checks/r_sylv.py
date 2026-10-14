# Referee: exact m_max at n = 2^k (k=12..16) and the Sylvester k-column cheapness window (Prop 6.5(c)).
import math
from mpmath import mp, loggamma, log, mpf
mp.dps = 60
def mmax(d):
    # locate via log-gamma then confirm exactly
    target = d*log(3) - log(8)
    lo, hi = 1, 2*d
    f = lambda m: loggamma(m+d) - loggamma(d+1) - loggamma(m) <= target
    while lo < hi:
        mid = (lo+hi+1)//2
        if f(mid): lo = mid
        else: hi = mid-1
    m = lo
    P = 3**d
    ok = lambda m: 8*math.comb(m+d-1, d) <= P
    assert ok(m) and not ok(m+1)
    return m
for k in range(12, 17):
    n = 2**k; d = n-1; r = n-1
    w = math.sqrt(8*d*math.log(4*r)); j = int(w//2)
    m = mmax(d)
    print(f"k={k} n={n} w={w:.3f} j={j} ceil(4w)={math.ceil(4*w)} m_max={m} (exact) 2jk={2*j*k} k-cols cheap inside window: {2*j*k < m}")
# k >= 23 analytic tail: 2jk <= k w and m_max > 0.6417 d - 2.2137
for k in range(17, 80):
    n = 2**k; d = n-1; w = math.sqrt(8*d*math.log(4*(n-1)))
    assert k*w < 0.6417*d - 2.2137, k
print("k*w < 0.6417 d - 2.2137 for 17 <= k < 80 (and k^{3/2} 2^{-k/2} decreasing beyond)")
