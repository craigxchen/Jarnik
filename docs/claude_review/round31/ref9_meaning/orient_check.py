"""Own check of the projection bound (IX.1) and of the orientation reading of (IX.2).
Exact integer arithmetic.  For every circle x^2+y^2=N, N<=NMAX, every window of consecutive
lattice points with chord constant max|zi-zj|/N^(1/4) <= CMAX, and every Gaussian b with
|b|^2 <= BMAX, compute X_b = diameter of Re(conj(b) z_i).  Report the largest values of
  r1 = M^(3/2)/(2+X_b)   and   r2 = M^2 exp(-sum_{p<=M, p|N|b|^2, p=1(4)} log p/p)/(2+X_b),
and compare with the trivial count bound M <= 2(X_b+1) and the record's round-4 Lemma 4.2(b)."""
import sys, math
from math import isqrt
NMAX=int(sys.argv[1]) if len(sys.argv)>1 else 200000
CMAX=float(sys.argv[2]) if len(sys.argv)>2 else 6.0
BMAX=50
pts={}
for x in range(-isqrt(NMAX),isqrt(NMAX)+1):
    for y in range(-isqrt(NMAX-x*x),isqrt(NMAX-x*x)+1):
        n=x*x+y*y
        if n: pts.setdefault(n,[]).append((x,y))
bs=[(a,b) for a in range(0,8) for b in range(-7,8) if 0<a*a+b*b<=BMAX and (a>0 or b>0)]
def primes(n):
    s=bytearray([1])*(n+1); s[0]=s[1]=0
    for i in range(2,isqrt(n)+1):
        if s[i]: s[i*i::i]=bytearray(len(s[i*i::i]))
    return [i for i in range(n+1) if s[i]]
PR=[p for p in primes(1000) if p%4==1]
best1=(0,None); best2=(0,None); cnt=0; maxM=0
for N,L in pts.items():
    if len(L)<6: continue
    L.sort(key=lambda z: math.atan2(z[1],z[0]))
    k=len(L); q=N**0.25
    for s in range(k):
        win=[L[s]]
        for t in range(1,k):
            z=L[(s+t)%k]
            if max(math.hypot(z[0]-w[0],z[1]-w[1]) for w in win)>CMAX*q: break
            win.append(z)
        M=len(win)
        if M<3: continue
        maxM=max(maxM,M); cnt+=1
        for (a,b) in bs:
            xs=[a*x+b*y for (x,y) in win]   # Re(conj(a+bi)(x+iy)) = ax+by
            X=max(xs)-min(xs)
            D=N*(a*a+b*b)
            loss=sum(math.log(p)/p for p in PR if p<=M and D%p==0)
            r1=M**1.5/(2+X); r2=M*M*math.exp(-loss)/(2+X)
            assert M<=2*(X+1)
            if r1>best1[0]: best1=(r1,(N,M,(a,b),X))
            if r2>best2[0]: best2=(r2,(N,M,(a,b),X))
print("NMAX",NMAX,"CMAX",CMAX,"windows with M>=3:",cnt,"max M",maxM)
print("max M^1.5/(2+X_b):",best1)
print("max M^2 exp(-loss)/(2+X_b):",best2)
