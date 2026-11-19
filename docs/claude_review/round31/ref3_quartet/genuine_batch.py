"""Many genuine same-unit quartets on short arcs: verify (III.1), the cofactor bound, and that the
genuine (x,z) lies in the congruence lattice A x + C z = 0 mod B (so the shortest-vector test is passed)."""
import math, random
from itertools import combinations
from math import gcd, prod
from genuine import points, analyse, det, gnorm
random.seed(31)
split=[p for p in range(5,400) if p%4==1 and all(p%q for q in range(2,int(p**.5)+1))]
tot=0; circles=0
while tot<3000 and circles<400:
    ps=random.sample(split,random.randint(5,8)); n=prod(ps)*random.choice([1,5,13,25])
    pts=points(n); N0=n
    # group by unit class relative to a fixed base orientation: use analyse on whole circle
    N,f,data=analyse(pts) if False else (None,None,None)
    from genuine import factorint, two_sq, vpi, gconj, gmul, gpow, UNITS
    f=factorint(n); pis={p:two_sq(p) for p in f}
    cls={}
    for z in pts:
        a={p:vpi(z,pis[p]) for p in f}
        base=(1,0)
        for p in f: base=gmul(base,gmul(gpow(pis[p],a[p]),gpow(gconj(pis[p]),f[p]-a[p])))
        eps=next(u for u in UNITS if gmul(u,base)==z)
        P=(1,0)
        for p in f: P=gmul(P,gpow(pis[p],a[p]))
        cls.setdefault(eps,[]).append((math.atan2(z[1],z[0]),z,a,P))
    circles+=1
    width=60*n**-0.25
    for eps,L in cls.items():
        L.sort()
        for s0 in range(len(L)):
          for q in ([L[s0]]+list(c) for c in combinations([u for u in L[s0+1:] if u[0]-L[s0][0]<=width],3)):
              q=sorted(q); span=q[-1][0]-q[0][0]
              if span>width or span<=0: continue
              # primitivity of the 4-tuple not required for identity; reduce by common exponent anyway
              ex=[r[2] for r in q]; Ps=[r[3] for r in q]
              g={};t={};d={}
              for i,j in combinations(range(4),2):
                  g[i,j]=prod(p**min(ex[i][p],ex[j][p]) for p in f)
                  D=det(Ps[i],Ps[j]); assert D%g[i,j]==0 and D!=0
                  t[i,j]=D//g[i,j]; d[i,j]=gnorm(Ps[i])*gnorm(Ps[j])//g[i,j]**2
                  assert abs(t[i,j])<=math.sqrt(d[i,j])*span/2*(1+1e-9)
              raw=[g[0,1]*g[2,3],g[0,2]*g[1,3],g[0,3]*g[1,2]]; G=gcd(*raw); A,B,C=[r//G for r in raw]
              x,y,z=t[0,1]*t[2,3],t[0,2]*t[1,3],t[0,3]*t[1,2]
              assert A*x-B*y+C*z==0 and x and y and z
              if B>1: assert (A*x+C*z)%B==0
              tot+=1
print('circles',circles,'same-unit quartets with span <= 60 n^(-1/4):',tot,'all pass (III.1), cofactor bound, congruence')
