# Recompute Prop 6.6 numbers and Sec 7 rank claims independently.
import sys; sys.path.insert(0,'.')
from gauss import *
from math import isqrt, atan2, sqrt
from fractions import Fraction
def hsum0(u): return sum(h_ratio(u[x],u[0]) for x in range(1,len(u)))
for (A,n,K) in [(10,20,2000),(2,30,3000),(3,8,50)]:
    g=(A,1); gp=[(1,0)]; 
    for _ in range(K): gp.append(mul(gp[-1],g))
    F1=[gp[k] for k in range(n)]+[gp[K]]
    F2=[gp[0]]+[gp[k] for k in range(K-n+1,K+1)]
    hl1,G1=h_line(F1); hl2,G2=h_line(F2)
    print(f"(A,n,K)=({A},{n},{K}): h(F1)/h(ell)={h_proj(F1)/hl1:.6f} [bound {K/(K+n-1):.6f}] G1={G1}; "
          f"hsum0(F2)/h(ell)={hsum0(F2)/hl2:.4f} [bound {n*(K-n+1)/(2*K):.4f}, n/2={n/2}] log2(G2)={G2.bit_length()-1}",flush=True)
# Section 7 data: exponent rank of best k-windows
def pts(N):
    P=[]
    for x in range(-isqrt(N),isqrt(N)+1):
        y2=N-x*x; y=isqrt(y2)
        if y*y==y2:
            P.append((x,y)); 
            if y: P.append((x,-y))
    return sorted(set(P),key=lambda p: atan2(p[1],p[0]))
def factor(n):
    f={};d=2
    while d*d<=n:
        while n%d==0: f[d]=f.get(d,0)+1; n//=d
        d+=1
    if n>1: f[n]=f.get(n,0)+1
    return f
def gprime(p):  # a+bi with a^2+b^2=p
    for a in range(1,isqrt(p)+1):
        b=isqrt(p-a*a)
        if a*a+b*b==p: return (a,b)
def vexp(z,pi):  # valuation of z at Gaussian prime pi
    e=0
    while True:
        t=mul(z,conj(pi)); n=norm(pi)
        if t[0]%n or t[1]%n: return e
        z=(t[0]//n,t[1]//n); e+=1
def rank(rows):
    M=[[Fraction(v) for v in r] for r in rows]; rk=0; cols=len(M[0]) if M else 0
    for c in range(cols):
        piv=next((i for i in range(rk,len(M)) if M[i][c]!=0),None)
        if piv is None: continue
        M[rk],M[piv]=M[piv],M[rk]
        for i in range(len(M)):
            if i!=rk and M[i][c]!=0:
                f=M[i][c]/M[rk][c]; M[i]=[a-f*b for a,b in zip(M[i],M[rk])]
        rk+=1
    return rk
for N,k in [(1176852625,10),(3159797225,7),(3159797225,6)]:
    P=pts(N); m=len(P); f=factor(N); sp=[p for p in f if p%4==1]
    best=None
    for s in range(m):
        w=[P[(s+j)%m] for j in range(k)]
        dx=w[-1][0]-w[0][0]; dy=w[-1][1]-w[0][1]; ch=sqrt(dx*dx+dy*dy)
        if best is None or ch<best[0]: best=(ch,w)
    ch,w=best; R=sqrt(N)
    z=w; g=z[0]
    for x in z[1:]: g=ggcd(g,x)
    pis=[gprime(p) for p in sp]
    rows=[[vexp(x,pi)-vexp(z[0],pi) for pi in pis] for x in z[1:]]
    print(f"N={N}={f} k={k}: best chord/sqrtR={ch/sqrt(R):.4f} gcd norm={norm(g)} split primes={len(sp)} rank B={rank(rows)} (M-1={k-1})",flush=True)
