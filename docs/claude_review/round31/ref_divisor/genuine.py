# Genuine clusters: normalize, profile, edge labels, floor reps, triangle mechanism, line reconstruction.
import sys, math
from math import isqrt, gcd
from itertools import combinations
from sympy import factorint
from fractions import Fraction
def mul(a,b): return (a[0]*b[0]-a[1]*b[1],a[0]*b[1]+a[1]*b[0])
def cj(a): return (a[0],-a[1])
def nrm(a): return a[0]**2+a[1]**2
def ggcd(z,w):
    while w!=(0,0):
        n=nrm(w); x,y=mul(z,cj(w)); q=((2*x+n)//(2*n),(2*y+n)//(2*n)); a,b=mul(q,w); z,w=w,(z[0]-a,z[1]-b)
    return z
def div(z,w):
    n=nrm(w); x,y=mul(z,cj(w)); assert x%n==0 and y%n==0; return (x//n,y//n)
def points(N):
    P=[]
    for x in range(isqrt(N)+1):
        y2=N-x*x; y=isqrt(y2)
        if y*y==y2:
            for sx in (1,-1):
                for sy in (1,-1): P.append((sx*x,sy*y))
    return sorted(set(P),key=lambda z:math.atan2(z[1],z[0]))
def best_cluster(N,k):
    P=points(N); n=len(P); best=None
    for i in range(n):
        Q=[P[(i+j)%n] for j in range(k)]
        a0=math.atan2(Q[0][1],Q[0][0]); a1=math.atan2(Q[-1][1],Q[-1][0])
        ang=(a1-a0)%(2*math.pi)
        if best is None or ang<best[0]: best=(ang,Q)
    ang,Q=best; return Q, ang*math.sqrt(N)/N**0.25
def analyse(N,k):
    Q,C=best_cluster(N,k)
    g=Q[0]
    for z in Q[1:]: g=ggcd(g,z)
    Z=[div(z,g) for z in Q]; Nn=nrm(Z[0])
    print(f"N={N} k={k} C={C:.6f}; normalized N'={Nn}={factorint(Nn)}")
    m=len(Z)
    lab={}; passes={}
    for i,j in combinations(range(m),2):
        gij=ggcd(Z[i],Z[j]); d=Nn//nrm(gij); lab[(i,j)]=d
        # genuine half-angle rep: z_i/z_j = A/conj A
        s=isqrt(d); r=d-s*s; t=isqrt(r)
        passes[(i,j)]= (r>0 and t*t==r and gcd(s,t)==1 and 256*t**4*Nn<=d*d)
    print("  distinct edge labels",len(set(lab.values())),"of",len(lab),"; edges passing primitive floor test",sum(passes.values()))
    # genuine triangle check: for each triangle, product square and the sign mechanism (Q_eps<=sqrt N for the u=0 pattern?)
    sq=0
    for i,j,l in combinations(range(m),3):
        P_=lab[(i,j)]*lab[(j,l)]*lab[(i,l)]; sq+= isqrt(P_)**2==P_
    print("  triangles with square label product",sq,"of",len(list(combinations(range(m),3))))
    # cliques of passing edges
    best=[]
    for r in range(m,1,-1):
        for S in combinations(range(m),r):
            if all(passes[(a,b)] for a,b in combinations(S,2)): best=S;break
        if best: break
    print("  largest all-pass sub-clique",best)
    return Z,Nn,lab,passes
for N,k in ((1176852625,10),(32988780417125,8),(235370525,7),(276108509,3)):
    analyse(N,k)
print("---- unit-class diagnosis")
def halfrep(zi,zj,N):
    re,im=mul(zi,cj(zj)); g=gcd(N+re,im); return (N+re)//g, im//g
for N,k in ((276108509,3),(1176852625,10),(32988780417125,8)):
    Q,C=best_cluster(N,k); g=Q[0]
    for z in Q[1:]: g=ggcd(g,z)
    Z=[div(z,g) for z in Q]; Nn=nrm(Z[0])
    for i,j in combinations(range(len(Z)),2):
        s,t=halfrep(Z[j],Z[i],Nn); d=Nn//nrm(ggcd(Z[i],Z[j]))
        if N==276108509 or (s*s+t*t==d and 256*t**4*Nn<=d*d):
            print(N,(i,j),"half-angle rep s,t=",s,t,"s^2+t^2 / d_ij =",Fraction(s*s+t*t,d))
