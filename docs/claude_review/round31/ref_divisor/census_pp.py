# parents with exponents in {0,1,2} on first n split primes; imprimitive candidates allowed.
import sys
from math import isqrt, gcd
from itertools import combinations, product
from collections import defaultdict
n=int(sys.argv[1]); E=2
P=[];q=5
def isp(x): return x>1 and all(x%k for k in range(2,isqrt(x)+1))
while len(P)<n:
    if q%4==1 and isp(q): P.append(q)
    q+=1
if len(sys.argv)>2: P=[int(x) for x in sys.argv[2].split(',')]; n=len(P)
cands=[]
for ev in product(range(E+1),repeat=n):
    d=1
    for p,e in zip(P,ev): d*=p**e
    if d==1: continue
    s=isqrt(d); r=d-s*s
    if r==0: continue
    t=isqrt(r)
    if t*t!=r or 256*t**4>d: continue
    cands.append((ev,d,s,t))
print("primes",P,"candidates",len(cands))
par=defaultdict(list)
for ev,d,s,t in cands:
    lim=d//(256*t**4)
    stack=[(0,1,list(ev))]
    while stack:
        k,m,cur=stack.pop()
        if k==n:
            par[tuple(cur)].append((ev,d,s,t)); continue
        f=1
        for extra in range(E-ev[k]+1):
            if m*f>lim: break
            c2=cur[:]; c2[k]=ev[k]+extra
            stack.append((k+1,m*f,c2))
            f*=P[k]
def mul(a,b): return (a[0]*b[0]-a[1]*b[1],a[0]*b[1]+a[1]*b[0])
cj=lambda a:(a[0],-a[1])
inc=0;nsqf=0;nimp=0;trip=set();bad=0
for Nv,v in par.items():
    for x,y,z in combinations(v,3):
        ok=True
        for k in range(n):
            a=sorted((x[0][k],y[0][k],z[0][k]))
            if a[0]+a[1]!=a[2] or a[2]>Nv[k]: ok=False;break
        if not ok: continue
        inc+=1
        if any(e==2 for e in Nv): nsqf+=1
        if any(gcd(w[2],w[3])>1 for w in (x,y,z)): nimp+=1
        key=tuple(sorted((x[1],y[1],z[1])));
        A=[(w[2],w[3]) for w in (x,y,z)]
        good=False
        for k in range(3):
            i,j=[l for l in range(3) if l!=k]
            pr=mul(mul(A[i],A[j]),cj(A[k]))
            if pr[1]==0 and pr[0]>0: good=True
        if not good: bad+=1
        trip.add(key)
print("incidences",inc,"at nonsquarefree parents",nsqf,"with imprimitive factor",nimp,"distinct triples",len(trip),"angle failures",bad)
