import sys, json, math
import numpy as np
from reflib import *
from slabbnb import bnb
name,K,q0,Fsrc=sys.argv[1],int(sys.argv[2]),int(sys.argv[3]),sys.argv[4]; lift=[int(x) for x in sys.argv[5].split(",")]
thr=json.load(open(Fsrc)); P=thr["P"]; r=len(P); F=sorted(thr["forced"]); pis=[gprime(p) for p in P]
m=m_of(q0); n=m//4; L,iota,_=genuine_L(P,q0); k=(iota//n)%4
V=monos(r,K,half=True); A=np.array(V,dtype=float); a=np.zeros((len(V),4))
for i,v in enumerate(V):
    Dv=Dvals(Zv(pis,v)); N=norm_v(P,v); d=sum(x*y for x,y in zip(v,lift))
    for s in range(4):
        g=1
        for q in F:
            if q!=q0:
                t=val(Dv[s],q)
                if t: g*=q**t
        if d%n==0 and ((d%m)//n)*(1 if k==1 else 3)%4==s: g*=q0**(1+val(d//n,q0))
        a[i,s]=(1 if s%2==0 else 2**-0.5)*g/math.sqrt(N)
st,phi,nodes=bnb(A,a,r,cap=400000)
print(name,K,q0,"slab BnB:",st,nodes, list(phi) if phi is not None else "")
