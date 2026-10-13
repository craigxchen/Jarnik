# Exact check of the order m in determinant.md Thm 2.1 against (n-1)T and w (README §8.3 remark).
from math import comb, log, sqrt
def m_of(n, T):
    bound = (2*T+1)**(n-1)  # compare 8*binom <= (2T+1)^(n-1)
    lo, hi = 0, 10*n*n
    while lo < hi:
        mid = (lo+hi+1)//2
        if 8*comb(mid+n-2, n-1) <= bound: lo = mid
        else: hi = mid-1
    return lo
for n in [64, 128, 500, 2137, 4278]:
    T = n; m = m_of(n, T)
    for r in [n]:
        w = T*sqrt(8*(n-1)*log(4*r))
        print(f"n={n} m={m} m/(n(n-1))={m/(n*(n-1)):.4f} (n-1)T/m={(n-1)*T/m:.3f} w={w:.0f} m/(2w)={m/(2*w):.3f} n*w/m={n*w/m:.1f}")
