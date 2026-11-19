# Independent exact check of Astra VI.4 four-point Pell families.
import sys, math
from fractions import Fraction
from itertools import combinations
from sympy import factorint
from mpmath import mp, mpf, atan2, log
mp.dps=1500
def mul(a,b): return (a[0]*b[0]-a[1]*b[1], a[0]*b[1]+a[1]*b[0])
def cj(a): return (a[0],-a[1])
def N(a): return a[0]*a[0]+a[1]*a[1]
def ggcd(a,b):
    while b!=(0,0):
        n=N(b); p=mul(a,cj(b))
        q=((2*p[0]+n)//(2*n),(2*p[1]+n)//(2*n))
        r=(a[0]-mul(q,b)[0],a[1]-mul(q,b)[1]); assert N(r)<n
        a,b=b,r
    return a
def gdiv(a,b):
    n=N(b); p=mul(a,cj(b)); assert p[0]%n==0 and p[1]%n==0; return (p[0]//n,p[1]//n)
def Apow(v,m):  # A=[[5,2],[2,1]] (det 1), negative m allowed
    x,y=v
    if m>=0:
        for _ in range(m): x,y=5*x+2*y,2*x+y
    else:
        for _ in range(-m): x,y=x-2*y,-2*x+5*y
    return (x,y)
def words(blocks, pre0=(1,0), pre=(0,-1)):
    pts=[]
    for fl in ((),(0,1),(0,2),(1,2)):
        p=pre if fl else pre0
        for j,z in enumerate(blocks): p=mul(p, z if j in fl else cj(z))
        pts.append(p)
    return pts
def analyse(pts):
    H=N(pts[0]); assert all(N(p)==H for p in pts); assert len(set(pts))==4
    g=pts[0]
    for p in pts[1:]: g=ggcd(g,p)
    q=[gdiv(p,g) for p in pts]; Hp=N(q[0])
    d2=max(N((a[0]-b[0],a[1]-b[1])) for a,b in combinations(q,2))
    # exponent log(Delta)/log(H), Delta = diam/sqrt(H) (primitive)
    e=(log(mpf(d2))/2 - log(mpf(Hp))/2)/log(mpf(Hp))
    # arc constant: angular spread * R' / R'^(1/3)
    angs=[atan2(mpf(z[1]),mpf(z[0])) for z in q]
    Rp=mp.sqrt(mpf(Hp)); spread=max(angs)-min(angs)
    K13=spread*Rp/Rp**(mpf(1)/3); K25=spread*Rp/Rp**(mpf(2)/5)
    return g,Hp,e,(K13,K25),q
def cutprofile(q):
    # log-weight of each cut (partition of the 4 marks) over primes of Hp; returns dict partition->log weight / log Hp
    Hp=N(q[0]); f=factorint(Hp); W={}
    for p,a in f.items():
        # find Gaussian prime pi over p
        if p==2:
            pis=[(1,1)]
        elif p%4==3: pis=[(p,0)]
        else:
            # find pi
            import sympy
            x=sympy.sqrt_mod(-1,p); pi=ggcd((p,0),(x,1)); pis=[pi,cj(pi)]
        for pi in pis:
            vals=[]
            for z in q:
                v=0;w=z
                while True:
                    n=N(pi); pr=mul(w,cj(pi))
                    if pr[0]%n or pr[1]%n: break
                    w=(pr[0]//n,pr[1]//n); v+=1
                vals.append(v)
            # layer decomposition: for t=1..max, set S_t={k: vals[k]>=t}
            for t in range(1,max(vals)+1):
                S=frozenset(k for k in range(4) if vals[k]>=t)
                if 0<len(S)<4:
                    key=min(S, frozenset(range(4))-S, key=lambda s:(len(s),sorted(s)))
                    W[key]=W.get(key,0)+math.log(N(pi))
    L=2*math.log(Hp)  # each block counted via pi and conj(pi)
    return {tuple(sorted(k)):round(v/L,4) for k,v in sorted(W.items(),key=lambda kv:-kv[1])}
if __name__=="__main__":
    for rates,seeds in [((1,1,1),((3,1),(4,1),(5,1))),((1,1,1),((2,1),(7,3),(1,3))),((6,7,7),((3,1),(4,1),(5,1))),((6,7,7),((2,1),(7,3),(1,3))),((6,7,7),((1,0),(1,2),(3,7)))]:
        print("rates",rates,"seeds",seeds)
        for n in (1,2,3,5,8,12,20):
            bl=[Apow(v,r*n) for v,r in zip(seeds,rates)]
            g,Hp,e,K,q=analyse(words(bl))
            line=f"  n={n:3d} N(gcd)={N(g):6d} bits(Hp)={Hp.bit_length():5d} logDelta/logH={float(e):.5f} arc/R'^(1/3)={float(K[0]):.4f} arc/R'^(2/5)={float(K[1]):.4f}"
            if Hp.bit_length()<330: line+=" cuts="+str(cutprofile(q))
            print(line)
