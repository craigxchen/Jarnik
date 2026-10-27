# Remark 4.4 recheck with own code; lambda = min atom = min(p/2, 1-p) (author uses p/2 throughout)
from mpmath import mp, mpf, log, binomial, linspace
mp.dps=40
q=3251; j=257; m=2057; P=(q+1)//2
def G(p,correct):
    lam = min(p/2,1-p) if correct else p/2
    tail=binomial(P,j+1)*(p/2)**(j+1)*(1-p)**(P-j-1)
    return m*log(3)+(m/mpf(2))*log(1/lam)+log(tail)/2
for correct in (False,True):
    best=None
    for e in range(0,2001):
        p=mpf(10)**(-mpf(e)/100)
        if p>=1: continue
        g=G(p,correct)
        if best is None or g<best[0]: best=(g,p)
    # refine near 1
    for t in range(1,2000):
        p=1-mpf(t)/2000
        g=G(p,correct)
        if g<best[0]: best=(g,p)
    print('lambda=min atom' if correct else 'lambda=p/2', 'min G =',float(best[0]),'at p=',float(best[1]))
