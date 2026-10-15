import sys, math
import numpy as np
from reflib import *
from rverify import R
from slab import setup
from slabbnb import bnb
fn=sys.argv[1]
W=R(open(fn,"rb")).load(); P,K=W["P"],W["K"]; r=len(P)
labs={q:[int(x) for x in W["labels"][q]["lifts"][-1]] for q in W["labels"]}
V,A,a,dem=setup(P,K,labs)
# genuine sanity on same moduli
pis=[gprime(p) for p in P]; ag=np.zeros_like(a)
for i,v in enumerate(V):
    Dv=Dvals(Zv(pis,v)); N=norm_v(P,v)
    for t in range(4):
        g=1
        for q in labs:
            e=val(Dv[t],q)
            if e: g*=q**e
        ag[i,t]=(1 if t%2==0 else 2**-0.5)*g/math.sqrt(N)
print("genuine demands, BnB:", bnb(A,ag*(1-1e-7),r,cap=200000)[0::2])
keep=list(range(len(V)))
st=bnb(A[keep],a[keep],r)[0]; print("full:",st)
for i in list(keep):
    k2=[j for j in keep if j!=i]
    if bnb(A[k2],a[k2],r,cap=200000)[0].startswith("INFEASIBLE"): keep=k2
print("minimal infeasible subsystem:", [(V[i], [int(round(x)) for x in (a[i]*np.sqrt(norm_v(P,V[i]))/np.array([1,2**-0.5,1,2**-0.5]))]) for i in keep])
