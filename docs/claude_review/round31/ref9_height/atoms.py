# Floating-point sanity test of the IX.4 atom bound sup_b P(sum w_i X_i=b) vs 1/(A*||w||_2) (and 1/(A sqrt(q h)))
import numpy as np, random, math
from math import gcd
from functools import reduce
def chi(x,q):
    x%=q
    return 0 if x==0 else (1 if pow(x,(q-1)//2,q)==1 else -1)
def paley(q): return np.array([[chi(j-i,q) for j in range(q)] for i in range(q)],dtype=np.int64)
def dg(A):
    R=int(12*A)+2; n=np.arange(-R,R+1); p=np.exp(-n**2/(2*A*A)); return p/p.sum()
def atom(w,A):
    k=dg(A); R=(len(k)-1)//2; pmf=np.array([1.0])
    for wi in w:
        if wi==0: continue
        ker=np.zeros(2*R*abs(wi)+1); ker[::abs(wi)]=k
        pmf=np.convolve(pmf,ker)
    return pmf.max()
def support_counts(w,h):
    out={}
    for n in [4]+[p for p in range(3,h+1) if all(p%r for r in range(2,p))]:
        out[n]=sum(1 for x in w if x%n)
    return out
random.seed(1)
for q in [103,199]:
    K=paley(q); rho=0.9
    Acls=[i for i in range(q) if K[i,0]==1]
    S=set(i for i in Acls if random.random()<rho); Sc=[i for i in range(q) if i not in S]
    labels=list(range(1,q))
    print(f'q={q} |S|={len(S)} |Sc|={len(Sc)}')
    for h in [2,4,6,8,12]:
        worst=0;rec=None
        for trial in range(6):
            while True:
                t=h//2; pos=[random.choice(labels) for _ in range(t)]; neg=[random.choice(labels) for _ in range(t)]
                v=np.zeros(q,dtype=np.int64)
                for j in pos: v[j]+=1
                for j in neg: v[j]-=1
                if np.abs(v).sum()==h and reduce(gcd,[int(x) for x in v if x])==1: break
            w=[int(x) for x in (K@v)[Sc]]
            g=reduce(gcd,[abs(x) for x in w if x])
            for A in [0.6,1.0,2.0,4.0]:
                at=atom(w,A); nw=math.sqrt(sum(x*x for x in w))
                r1=at*A*nw*math.sqrt(2*math.pi); r2=at*A*math.sqrt(q*h)
                if r1>worst: worst=r1; rec=(A,g,r2,support_counts(w,h))
        print(f'  h={h}: max over samples of atom*A*||w||*sqrt(2pi) = {worst:.3f}  (A,gcd,atom*A*sqrt(qh),supp-counts)={rec}')
# adversarial w (not from Paley): all multiples of 3 except n1 entries equal 1
print('adversarial: w = 3 on N entries, 1 on n1 entries; ratio atom*A*||w||*sqrt(2pi)')
for n1 in [0,1,2,4,8,16]:
    w=[3]*60+[1]*n1
    for A in [0.6,1.0]:
        at=atom(w,A); nw=math.sqrt(sum(x*x for x in w))
        print(f'  n1={n1} A={A}: ratio={at*A*nw*math.sqrt(2*math.pi):.3f}  bound-shape A*Z/sqrt(n1)-type term')
