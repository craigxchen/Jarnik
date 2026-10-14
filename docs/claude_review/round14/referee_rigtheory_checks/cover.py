import itertools
from glib import *
from forced import forced_set
ps=[5,13,17,29]; r=4
def ball(h): return [v for v in itertools.product(range(-h,h+1),repeat=r) if sum(map(abs,v))<=h]
for K in (3,4,5,6):
    h=(K-1)//2; F,_=forced_set(ps,K,1200); Ds_=[]
    if (K-1)%2==0: Ds_=[ball(h)]
    else:
        for j in range(r):
            for sg in (1,-1):
                e=tuple(sg if t==j else 0 for t in range(r)); B=ball(h)
                Ds_.append(list(set(B)|{tuple(a+b for a,b in zip(v,e)) for v in B}))
    cov=[]
    for q in F:
        m,L=genuine(ps,q); n=m//4
        if any(len({sum(a*b for a,b in zip(v,L))%n for v in D})==n for D in Ds_): cov.append(q)
    print("K",K,"|D|",len(Ds_[0]),"covered forced moduli:",cov)
