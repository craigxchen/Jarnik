import math
from math import comb, lgamma
from scans import params
def m_max(n):
    d = n-1; L3 = d*math.log(3)
    # float estimate then exact refine
    f = lambda m: math.log(8) + lgamma(m+d) - lgamma(d+1) - lgamma(m)
    lo, hi = 1, d
    while lo < hi:
        mid = (lo+hi+1)//2
        if f(mid) <= L3: lo = mid
        else: hi = mid-1
    m = lo; T = 3**d
    while 8*comb(m+d-1, d) > T: m -= 1
    while 8*comb(m+d, d) <= T: m += 1
    return m
if __name__ == '__main__':
  for rk in ['n','nlnn','n2']:
    fails = [n for n in range(2000, 6500) if params(n, rk)[3] > m_max(n)]
    print(rk, 'last n in [2000,6500) with ceil(4w) > m_max:', fails[-1] if fails else None, flush=True)
