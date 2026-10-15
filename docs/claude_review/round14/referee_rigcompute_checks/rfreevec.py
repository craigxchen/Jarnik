# verify the author's stored FREE vectors and coverage of [3,Qscan] for a thr json
import sys, json, numpy as np
from reflib import *
for fn in sys.argv[1:]:
    d=json.load(open(fn)); P=d['P']; K=d['K']; r=len(P)
    V=np.array(monos(r,K,half=True),dtype=np.int64)
    forced=set(d['forced']); ab=set(d['aborted']); free={int(k):v for k,v in d['free'].items()}
    qs=[q for q in primerange(3,d['Qscan']+1) if q not in P]
    missing=[q for q in qs if q not in forced and q not in ab and q not in free]
    bad=[]
    for q,c in free.items():
        n=m_of(q)//4; nu=nu_par(P,q)
        if n%2==0 and any(x%2!=p for x,p in zip(c,nu)): bad.append((q,'par')); continue
        if np.any((V@np.array(c,dtype=np.int64))%n==0): bad.append((q,'coll'))
    print(fn.split('/')[-1],"K",K,"Qscan",d['Qscan'],"#forced",len(forced),"#aborted",len(ab),"#free",len(free),"missing",missing[:5],len(missing),"bad free",bad[:5],flush=True)
