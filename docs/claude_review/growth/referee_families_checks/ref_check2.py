import sys, itertools, math
sys.path.insert(0,'.')
import ref_check as rc   # reruns golden part on import (quick)
from decimal import Decimal as D, getcontext
getcontext().prec = 1000
mul,cj=rc.mul,rc.cj
# notes' translated Pell eight-point family, n = 1 mod 60
def pell(r):
    x,y=1,0
    for _ in range(r): x,y=9*x+20*y,4*x+9*y
    return x,y
def H(r):
    x,y=pell(r); return (x+y,2*y)
for n in (61,):
    bl=[H(n-2),H(n),H(n+2),H(n+4)]; pts=[]
    for w in itertools.product([0,1],repeat=4):
        p=sum(w)
        if p%2==0: continue
        z=rc.F if p==1 else cj(rc.F)
        for b,s in zip(bl,w): z=mul(z,b if s else cj(b))
        pts.append(z)
    C,N,k,gn,Ps=rc.cluster_C(pts)
    print("Pell8 n=%d k=%d gcdnorm=%d C=%s log10N=%.1f"%(n,k,gn,str(C)[:16],math.log10(N)))
    # is H_r associate to the golden J_{6r}?  compare norms with fib
    r=5; print("Norm H_5 =",rc.nrm(H(5))," f_61=",rc.fib[61], " Norm J_30=",rc.nrm(rc.J(30)))
