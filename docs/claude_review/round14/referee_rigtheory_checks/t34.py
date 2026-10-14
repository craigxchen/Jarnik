import random, math, itertools
from glib import *
ps=[5,13,17,29]; r=4
def ball(h):
    return [v for v in itertools.product(range(-h,h+1),repeat=r) if sum(map(abs,v))<=h]
random.seed(1)
for q,h,K in [(157,2,5),(41,1,3),(263,2,6),(79,1,4)]:
    m,L=genuine(ps,q); n=m//4; c=[x%n for x in L]
    D=ball(h)
    if K%2==0:  # double ball for even K (diameter 2h+1 = K-1)
        D=list(set(D)|set(tuple(a+(1 if j==0 else 0) for j,a in enumerate(v)) for v in D))
    cov=len({sum(a*b for a,b in zip(v,c))%n for v in D})==n
    DD={tuple(a-b for a,b in zip(x,y)) for x in D for y in D}
    S1=[tuple(1 if j==i else 0 for j in range(r)) for i in range(r)]; S1+= [tuple(-x for x in e) for e in S1]+[(0,)*r]
    DSD={tuple(a+e-b for a,e,b in zip(x,ee,y)) for x in D for y in D for ee in S1}
    W=[w for w in DD|DSD if any(w)]
    cls={tuple((u*x)%n for x in c) for u in range(1,n) if math.gcd(u,n)==1}
    fails=0; tested=0
    for _ in range(3000):
        cp=tuple(random.randrange(n) for _ in range(r))
        if cp in cls: continue
        tested+=1
        if not any(sum(a*b for a,b in zip(w,cp))%n==0 and sum(a*b for a,b in zip(w,c))%n!=0 for w in W): fails+=1
    print(q,"n",n,"|D|",len(D),"covered",cov,"tested",tested,"counterexamples",fails)
