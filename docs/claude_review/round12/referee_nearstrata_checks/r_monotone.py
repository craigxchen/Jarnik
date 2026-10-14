# Does the Cor 3.7 / n_B condition (proved orders) at fixed n hold for every r <= n once it holds at r = n?
import math, numpy as np
def Bval(n, r):
    d = n-1; w = math.sqrt(8*d*math.log(4*r)); j = math.floor(w/2); mp_ = math.ceil(4*w)
    if w < 10: return None
    N2 = n*(n-1)/2; b = (mp_-1)//(2*j); worst = -1e300
    for s in range(b+1):
        emax = (mp_-1)//2 - j*s
        if emax < 0: continue
        e = np.arange(emax+1, dtype=float)
        Q = 2**(s+1)*N2/((n-e)*(n-e-2**s))
        worst = max(worst, float(np.max((mp_-2*j*s-2*e)*Q)) - d)
    return worst
for n in [5185, 5500, 6000, 8000, 12000, 20000]:
    bad = [r for r in range(1, n+1) if (Bval(n, r) is not None and Bval(n, r) >= 0)]
    print(n, 'B(r=n)=', round(Bval(n,n),2), ' #r<=n failing:', len(bad), ' examples:', bad[:5], bad[-3:], flush=True)
