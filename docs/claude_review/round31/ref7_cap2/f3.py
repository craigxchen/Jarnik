# F_3[t] arithmetic with int arrays (index=degree)
import numpy as np
def tr(a):
    a=np.mod(np.asarray(a,dtype=np.int64),3); nz=np.nonzero(a)[0]
    return a[:nz[-1]+1] if len(nz) else np.zeros(1,np.int64)
def isz(a): return len(a)==1 and a[0]==0
def pmod(a,b):
    a=tr(a).copy(); b=tr(b); db=len(b)-1; inv=1 if b[-1]==1 else 2
    while len(a)-1>=db and not isz(a):
        d=len(a)-1; c=(a[-1]*inv)%3
        a[d-db:d+1]=(a[d-db:d+1]-c*b)%3
        a=tr(a)
    return a
def pdivmod(a,b):
    a=tr(a).copy(); b=tr(b); db=len(b)-1; inv=1 if b[-1]==1 else 2
    qd=max(len(a)-1-db,0); q=np.zeros(qd+1,np.int64)
    while len(a)-1>=db and not isz(a):
        d=len(a)-1; c=(a[-1]*inv)%3; q[d-db]=c
        a[d-db:d+1]=(a[d-db:d+1]-c*b)%3; a=tr(a)
    return tr(q),a
def gcd(a,b):
    a,b=tr(a),tr(b)
    while not isz(b): a,b=b,pmod(a,b)
    inv=1 if a[-1]==1 else 2
    return tr(a*inv)
def deriv(a):
    a=tr(a); return tr(a[1:]*np.arange(1,len(a))) if len(a)>1 else np.zeros(1,np.int64)
