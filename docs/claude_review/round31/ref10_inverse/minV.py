# Referee's own code: minimum owner-row weights over all integer V with D V = (H-D)/2
import sys, math, random
import numpy as np
from scipy.optimize import milp, LinearConstraint, Bounds
from sympy import Matrix
from sympy.matrices.normalforms import smith_normal_form
from sympy import ZZ

def chi(a,q):
    a%=q
    if a==0: return 0
    return 1 if pow(a,(q-1)//2,q)==1 else -1

def setup(q, S=None):
    I=list(range(1,q))
    if S is None: S=[i for i in I if chi(-i,q)==1]
    i0=S[0]
    a={(i,k):1+chi(k-i,q) for i in I for k in I}
    rows=[i for i in S if i!=i0]
    D=np.array([[a[(i,k)]-a[(i0,k)] for k in I] for i in rows],dtype=np.int64)
    Hm=np.array([[(k==i)-(k==i0) for k in I] for i in rows],dtype=np.int64)
    assert ((Hm-D)%2==0).all()
    Bm=(Hm-D)//2
    return I,S,i0,rows,D,Hm,Bm

def smith_ok(D):
    M=Matrix(D.tolist())
    sn=smith_normal_form(M,domain=ZZ)
    d=[sn[i,i] for i in range(min(sn.shape))]
    return all(abs(x)==1 for x in d), sorted(set(abs(x) for x in d))

TL=60
def solve_col(D,b,I,S,wts,owner_only=True,bnd=60):
    n=len(I); idx={k:t for t,k in enumerate(I)}
    # vars: v (n), t (n) with t>=|v|
    c=np.zeros(2*n)
    for t,k in enumerate(I):
        if k in S: c[n+t]=wts[t]
        elif not owner_only: c[n+t]=wts[t]*1.0
        else: c[n+t]=1e-4  # tiny tie-break to keep non-owners small
    A1=np.hstack([D,np.zeros_like(D)])
    E=np.eye(n)
    A2=np.hstack([E,-E]); A3=np.hstack([-E,-E])
    cons=[LinearConstraint(A1,b,b),LinearConstraint(A2,-np.inf,0),LinearConstraint(A3,-np.inf,0)]
    integ=np.concatenate([np.ones(n),np.zeros(n)])
    lb=np.concatenate([-bnd*np.ones(n),np.zeros(n)]); ub=np.concatenate([bnd*np.ones(n),bnd*np.ones(n)])
    r=milp(c,constraints=cons,integrality=integ,bounds=Bounds(lb,ub),options={'time_limit':TL})
    assert r.x is not None, r.message
    v=np.round(r.x[:n]).astype(int)
    assert (D@v==b).all()
    return v

def run(q,S=None,wts=None,owner_only=True):
    I,S,i0,rows,D,Hm,Bm=setup(q,S)
    n=len(I); m=len(S)
    if wts is None: wts=np.ones(n)
    ok,inv=smith_ok(D)
    V=np.zeros((n,n),dtype=int)  # V[j,k]
    for t,k in enumerate(I):
        V[:,t]=solve_col(D,Bm[:,t],I,S,wts,owner_only)
    # check D F = H with F = I+2V
    F=np.eye(n,dtype=int)+2*V
    assert (D@F==Hm).all()
    Sidx=[I.index(i) for i in S]
    roww=[int(np.dot(np.abs(V[s]),wts)) if np.allclose(wts,1) else float(np.dot(np.abs(V[s]),wts)) for s in Sidx]
    # parity bound: per column, min_c #{i in S: (c+B_ik) odd}, B row of i0 = 0
    Bfull={i0:np.zeros(n,dtype=int)}
    for r_,i in enumerate(rows): Bfull[i]=Bm[r_]
    par=0.0
    for t in range(n):
        cnt=[sum(wts[t] for i in S if (c+Bfull[i][t])%2==1) for c in (0,1)]
        par+=min(cnt)
    T=wts.sum(); Bq=(wts**2).sum(); Pm=wts.max()
    cs=0.5*(m*T-math.sqrt(m*q*Bq)-m*Pm)
    tot=sum(roww)
    # rows of all labels too
    allw=[int(np.abs(V[j]).sum()) for j in range(n)]
    return dict(q=q,m=m,smith_all1=ok,inv=inv,total_owner=tot,parity_bound=par,CS_bound=cs,
                per_row_mean=tot/m,half_T=T/2,min_row=min(roww),max_row=max(roww),
                rows=sorted(roww),maxentry=int(np.abs(V).max()),nonowner_total=sum(allw)-sum(int(np.abs(V[s]).sum()) for s in Sidx))

if __name__=="__main__":
    for q in [int(x) for x in sys.argv[1:]]:
        R=run(q)
        print(R,flush=True)
