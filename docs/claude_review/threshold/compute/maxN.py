from sizes2 import *
import sys
for M in range(3,10):
    for t in range(1,25):
        res=[]
        for s in range(t%2, min(t,4)+1, 2):
            p=(t+s)//2; n=(t-s)//2
            O=orbits(M,p,n); mus=[mu_of(v) for v in O]
            best=max((sum(kostka(lam,m) for m in mus),lam) for lam in partitions(M))
            res.append((s,best[0]))
        mx=max(r[1] for r in res)
        print(M,t,res,flush=True)
        if mx>40000: break
