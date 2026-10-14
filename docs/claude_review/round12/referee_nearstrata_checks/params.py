# Independent recomputation of regime parameters (exact where claimed exact).
import math, sys
from fractions import Fraction
from math import comb, isqrt
import mpmath as mp
mp.mp.dps = 50

def w_of(n, r):
    d = n-1
    return mp.sqrt(8*d*mp.log(4*r))

def j_mprime(n, r):
    w = w_of(n, r)
    j = int(mp.floor(w/2))
    mp_ = int(mp.ceil(4*w))
    return w, j, mp_

def m_max(n):
    # largest m with 8*binom(m+d-1,d) <= 3^d ; binom(m+d-1,d)=binom(m+d-1,m-1)
    d = n-1
    T = 3**d
    # binary search on m in [1, d]
    lo, hi = 1, d
    if 8*comb(lo+d-1, d) > T: return 0
    while lo < hi:
        mid = (lo+hi+1)//2
        if 8*comb(mid+d-1, d) <= T: lo = mid
        else: hi = mid-1
    return lo

if __name__ == "__main__":
    out = []
    for k in [14,15,16,17]:
        n = 2**k; r = n-1
        w, j, mp_ = j_mprime(n, r)
        mm = m_max(n)
        out.append(f"k={k} n={n} r={r} w={mp.nstr(w,12)} j={j} ceil4w={mp_} 2jk={2*j*k} m_max={mm} window_nonempty={2*j*k < mm}")
    print("\n".join(out))
