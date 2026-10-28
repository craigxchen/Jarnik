import sys; sys.path.insert(0,'.')
from gi import *
from fractions import Fraction
from itertools import combinations
import math, random
roots={'123':(0,1),'12':(-6,-4),'13':(-6,-3),'23':(-6,5),'1':(-7,6),'2':(-3,-2),'3':(-4,-3)}
def H(T,s): r=roots[T]; return (s-r[0],-r[1])
def pts(s):
    z0=(1,0)
    for T in roots: z0=gmul(z0,gconj(H(T,s)))
    Z=[z0]
    for x in '123':
        z=(1,0)
        for T in roots: z=gmul(z, H(T,s) if x in T else gconj(H(T,s)))
        Z.append(z)
    return Z
cuts=[frozenset(c) for c in [{0},{1},{2},{3},{0,1},{0,2},{0,3}]]
def analyse(s):
    Z=pts(s); g=ggcdlist(Z); zeta=[gexact(z,g) for z in Z]
    N=gnorm(zeta[0]); assert all(gnorm(z)==N for z in zeta)
    assert gnorm(ggcdlist(zeta))==1
    diam2=max(gnorm((a[0]-b[0],a[1]-b[1])) for a,b in combinations(zeta,2))
    logR=0.5*math.log(N)
    ratio=math.exp(0.5*math.log(diam2)-(3/7)*logR)
    # pair data
    d={};t={}
    for a,b in combinations(range(4),2):
        gab=ggcd(zeta[a],zeta[b]); diff=(zeta[a][0]-zeta[b][0],zeta[a][1]-zeta[b][1])
        c=gexact(diff,gab); n4=gnorm(c); assert n4%4==0; tt=isqrt(n4//4); assert tt*tt*4==n4
        t[(a,b)]=tt; d[(a,b)]=math.log(N)-math.log(gnorm(gab))
    # solve for log cut norms: d_ab = sum_{cut sep a,b} x_c ; sum x_c = W'
    import numpy as np
    A=[];y=[]
    for (a,b),v in d.items():
        A.append([1.0 if ((a in c)!=(b in c)) else 0.0 for c in cuts]); y.append(v)
    A.append([1.0]*7); y.append(math.log(N))
    x=np.linalg.solve(np.array(A),np.array(y))
    w=math.log(N)/7
    eta=max(abs(xi-w)/w for xi in x)
    return gnorm(g),ratio,max(t.values()),eta,[round(math.exp(v)) for v in x]
out=open('k4_lattice.out','w')
mx=[0,0,1e9,0,0]
random.seed(1)
S=[10,100,1000,10**4,10**6,10**9,10**12,10**15,10**20,10**25,10**30]+[random.randint(10**5,10**30) for _ in range(300)]
for s in S:
    gn,ratio,tmax,eta,cn=analyse(s)
    mx[0]=max(mx[0],gn); mx[1]=max(mx[1],ratio); mx[2]=min(mx[2],ratio); mx[3]=max(mx[3],tmax)
    if s in (10**6,10**15,10**30): print(s,gn,round(ratio,3),tmax,round(eta,5),file=out)
print('max gcd norm',mx[0],' ratio range',round(mx[2],3),round(mx[1],3),' max |t|',mx[3],file=out)
out.close(); print(open('k4_lattice.out').read())
