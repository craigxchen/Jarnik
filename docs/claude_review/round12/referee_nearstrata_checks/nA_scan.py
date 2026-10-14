import math
from fractions import Fraction as F
from scans import params
def A_val(n, rk, exact=False):
    d,w,j,mp_ = params(n, rk); N2 = n*(n-1)//2
    best = None
    for s in range(0,(mp_-1)//(2*j)+1):
        e = (mp_-1)//2 - j*s
        if e < 0: continue
        if exact: v = F(mp_*2**(s+1)*N2, (n-e)*(n-e-2**s)) - d
        else: v = mp_*2**(s+1)*N2/((n-e)*(n-e-2**s)) - d
        best = v if best is None or v > best else best
    return best
for rk, lo, hi in [('n', 300000, 4000000), ('nlnn', 400000, 2000000), ('n2', 700000, 2000000)]:
    last = None
    for n in range(lo, hi+1):
        if A_val(n, rk) > 0: last = n
    ex = [(m, float(A_val(m, rk, True))) for m in (last, last+1)] if last else None
    print(rk, 'last failing order-free (A>0) in', (lo,hi), ':', last, 'exact at last,last+1:', ex, flush=True)
