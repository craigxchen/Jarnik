# exact integer generators for the explicit types, from closed flats with 1-dim complement
import numpy as np, pickle
from flint import fmpq_mat, fmpz
from math import gcd
from fractions import Fraction
B=np.array([list(map(int,l.split())) for l in open('hadB.txt')])
want={(2,6),(4,6),(4,9),(6,2),(6,4),(7,5),(12,1)}
pools={k:set() for k in want}
for ln in open('r_flats_out.txt'):
    a=list(map(int,ln.split())); R,Z,s,t,dn=a[:5]
    if (s,t) not in want: continue
    assert dn==1
    cols=[x for x in range(12) if R>>x&1]
    rows=[B[0]]+[B[c+1] for c in range(11) if Z>>c&1]
    Mt=fmpq_mat([[int(r[x]) for x in cols] for r in rows])
    # nullspace via rref
    rr,rank=Mt.rref()
    piv=[];
    for i in range(rank):
        for j in range(len(cols)):
            if rr[i,j]!=0: piv.append(j); break
    free=[j for j in range(len(cols)) if j not in piv]; assert len(free)==1
    fc=free[0]; v=[Fraction(0)]*len(cols); v[fc]=Fraction(1)
    for i,pj in enumerate(piv): v[pj]=-Fraction(int(rr[i,fc].p),int(rr[i,fc].q))
    den=1
    for x in v: den=den*x.denominator//gcd(den,x.denominator)
    w=[int(x*den) for x in v]; g=0
    for x in w: g=gcd(g,abs(x))
    w=[x//g for x in w]
    z=[0]*12
    for j,x in zip(cols,w): z[x if False else j]=x
    nz=[x for x in z if x][0]
    if nz<0: z=[-x for x in z]
    z=tuple(z)
    # exact type check
    zz=np.array(z); assert zz.sum()==0
    assert np.count_nonzero(zz)==s
    tt=sum(1 for c in range(11) if (B[c+1]*zz).sum()!=0); assert tt==t,(z,tt,t)
    pools[(s,t)].add(z)
for k in sorted(pools): print(k,len(pools[k]),'max|entry|',max(max(abs(x) for x in z) for z in pools[k]))
pickle.dump({k:sorted(v) for k,v in pools.items()},open('r_pools.pkl','wb'))
