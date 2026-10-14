# independent check of Thm 3.4(i) (orbit count) and Thm 3.3 on several binary codes, all N_q = 2j+1
import itertools, sys
from rlib import *
def span(gens,b):
    C={tuple([0]*b)}
    for g in gens:
        C|={tuple((a+x)%2 for a,x in zip(c,g)) for c in C}
    return sorted(C)
def orbit_H(C,b,j,Dmax):
    E=list(itertools.product(range(2*j+1),repeat=b)); seen=set(); H=[0]*(2*j*b+1)
    for e in E:
        if e in seen: continue
        orb={}
        for s in C:
            f=tuple(2*j-x if si else x for x,si in zip(e,s))
            orb[f]=1
        for f in orb: seen.add(f)
        killed=False
        for s in C:
            if any(s) and all((2*j-x if si else x)==x for x,si in zip(e,s)):
                if (-1)**(j*sum(s))==-1: killed=True
        if not killed: H[min(sum(f) for f in orb)]+=1
    return H
cases=[((1,1,1),),((1,1,1,1),),((1,1,1,0,0),(0,0,1,1,1)),((1,1,0,1,0),(0,1,1,0,1)),((1,1,1,0),(0,1,1,1)),((1,0,0,1,1),(0,1,0,1,0),(0,0,1,0,1)),((1,1,0),),((1,0,1,1),)]
for gens in cases:
    b=len(gens[0])
    for j in (1,2):
        if (2*j+1)**b>4000: continue
        C=span(gens,b)
        G=list(itertools.product(range(-j,j+1),repeat=b))
        P=[y for y in G if all(sum(yy for yy,ss in zip(y,s) if ss)%2==0 for s in C)]
        cap=[2*j+1]*b
        hs=[];e=0
        while True:
            hs.append(h_func_exact(P,e,cap))
            if hs[-1]==len(P): break
            e+=1
        H=[hs[0]]+[hs[i]-hs[i-1] for i in range(1,len(hs))]
        Ho=orbit_H(C,b,j,0); 
        while Ho and Ho[-1]==0: Ho.pop()
        print(gens,"j",j,"|P|",len(P),"H_exact",H,"orbit",Ho,"match",H==Ho); sys.stdout.flush()
