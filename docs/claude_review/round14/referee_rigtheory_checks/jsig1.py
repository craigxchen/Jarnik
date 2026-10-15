"""Level-1 signatures (fast, chunked, numpy-unique) for moduli with lev=1. Same semantics as jsig.sigs(lev=1)."""
import math, itertools, numpy as np
from glib import *
from jsig import U1_orbit_reps, minimal, canon
def sigs1(ps, K, q):
    r=len(ps); pis=[two_sq(p) for p in ps]; V=[v for v in VK(r,K) if canon(v)]
    nz=np.array([Zv(pis,v)[0]**2+Zv(pis,v)[1]**2 for v in V],dtype=object)
    aff=np.array([[q*q <= int(nz[i])*(2 if s%2 else 1) for s in range(4)] for i in range(len(V))])
    m=q-1 if q%4==1 else q+1; n=m//4; nuv=nu(ps,q); Va=np.array(V,np.int64)
    m0,L=genuine(ps,q); cls=set()
    for eps in itertools.product([1,-1],repeat=r):
        for u in range(1,n):
            if math.gcd(u,n)==1: cls.add(tuple((u*e*x)%n for e,x in zip(eps,L)))
    ch=[np.arange(nuv[j],m,2) for j in range(r)]
    rest=np.array(np.meshgrid(*ch[1:],indexing='ij')).reshape(r-1,-1).T
    keys={}
    for a0 in U1_orbit_reps(m,nuv[0]):
        for k0 in range(0,len(rest),100000):
            blk=rest[k0:k0+100000]; grid=np.hstack([np.full((len(blk),1),a0),blk])
            t=(grid@Va.T)%m; coll=(t%n==0); s=(t//n)%4
            ok=((~coll)|aff[np.arange(len(V))[None,:],s]).all(axis=1)
            g=grid[ok]; key=np.where(coll[ok],s[ok]+1,0).astype(np.uint8)
            c=g%n; base=n**np.arange(r)[::-1]; code=c@base
            clsc=np.array(sorted(sum(x*b for x,b in zip(cc,base)) for cc in cls)) if n>1 else None
            inc=np.isin(code,clsc) if n>1 else np.ones(len(g),bool)
            for flag,mask in ((1,inc),(2,~inc)):
                if mask.any():
                    for row in np.unique(key[mask],axis=0):
                        sig=tuple((int(i)*4+int(row[i])-1,1) for i in np.nonzero(row)[0])
                        keys[sig]=keys.get(sig,0)|flag
    return keys
