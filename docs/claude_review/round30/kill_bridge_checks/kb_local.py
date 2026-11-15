# Independent local maxima of |g(t)|_l by brute-force l-adic sampling (no nearest-candidate argument).
import sys, itertools, sympy as sp
from fractions import Fraction as Fr
from kb_gauss import *
def vl(x,l):
    x=Fr(x); 
    if x==0: return 10**9
    v=0; n,d=x.numerator,x.denominator
    while n%l==0: n//=l; v+=1
    while d%l==0: d//=l; v-=1
    return v
def normg(rts,t0,masks=MASKS):
    g,den,ints=frac_gcd(values(rts,t0,masks))
    return Fr(g[0]**2+g[1]**2,den*den)
def badprimes(rts):
    P=set()
    for a,b in rts:
        for x in (a,b):
            P|=set(sp.factorint(Fr(x).denominator))
    allr=rts+[(a,-b) for a,b in rts]
    for (a,b),(c,d) in itertools.combinations(allr,2):
        n=(a-c)**2+(b-d)**2
        P|=set(sp.factorint(n.numerator)); P|=set(sp.factorint(n.denominator))
    P.discard(1); return sorted(P)
def localmax(rts,l,jmax,k,masks=MASKS):
    best=None;arg=None
    for j in range(0,jmax+1):
        M=l**(j+k)
        for a in range(M):
            t0=Fr(a,l**j)
            if t0.denominator!=l**j and j>0: continue
            v=vl(normg(rts,t0,masks),l)
            if best is None or v>best: best,arg=v,t0
    return best,arg
if __name__=="__main__":
    y=[Fr(-9,4),Fr(-2,3),Fr(1,4),Fr(8,3)]
    rts=roots_cube(y)
    P=badprimes(rts); print("P=",P)
    tot=Fr(1); import math
    res={}
    for l in P:
        if l==2: jm,k=6,9
        elif l==3: jm,k=4,6
        elif l<20: jm,k=1,3
        else: jm,k=0,2
        v,arg=localmax(rts,l,jm,k)
        res[l]=(v,arg); print(l, "max v_l(Ng)=",v, "|g|_l=",l**(v/2), "at",arg, flush=True)
    Gam=153335/1296
    logprod=sum(v*math.log(l)/2 for l,(v,a) in res.items())
    print("C_inf =",Gam/math.exp(logprod/2))
