# Attain the parity bound: owner entries forced to {-1,0,1} with the optimal parity pattern; feasibility => optimality
import sys, numpy as np
from minV import setup
from scipy.optimize import milp, LinearConstraint, Bounds
def col(D,b,I,S,par,TL):
    n=len(I)
    # v = 2y + par on owners (par in {0,1}), owner v in [-1,1]; non-owner integer
    # variables: v(n), y(n) integer; constraints D v = b; v_s - 2 y_s = par_s for owners
    N=2*n; c=np.zeros(N)
    A=[np.hstack([D,np.zeros_like(D)])]; lo=[b]; hi=[b]
    for s_,k in enumerate(I):
        if k in S:
            row=np.zeros(N); row[s_]=1; row[n+s_]=-2
            A.append(row[None,:]); lo.append([par[k]]); hi.append([par[k]])
    A=np.vstack(A); lo=np.concatenate(lo); hi=np.concatenate(hi)
    lb=np.full(N,-60.); ub=np.full(N,60.)
    for s_,k in enumerate(I):
        if k in S: lb[s_]=-1; ub[s_]=1
    r=milp(c,constraints=[LinearConstraint(A,lo,hi)],integrality=np.ones(N),bounds=Bounds(lb,ub),options={'time_limit':TL})
    if r.x is None: return None
    v=np.round(r.x[:n]).astype(int); assert (D@v==b).all(); return v
for q in map(int,sys.argv[1:]):
    I,S,i0,rows,D,Hm,Bm=setup(q); n=len(I)
    Bf={i0:np.zeros(n,dtype=int)}
    for r_,i in enumerate(rows): Bf[i]=Bm[r_]
    V=np.zeros((n,n),dtype=int); par_total=0; ok=True
    for t in range(n):
        opts=[]
        for c in (0,1):
            par={i:(c+Bf[i][t])%2 for i in S}; opts.append((sum(par.values()),par))
        opts.sort(key=lambda x:x[0])
        v=None
        for cnt,par in opts:
            if cnt>opts[0][0]: break
            v=col(D,Bm[:,t],I,S,par,120)
            if v is not None: break
        if v is None: ok=False; print('col',t,'not attained',flush=True); continue
        V[:,t]=v; par_total+=opts[0][0]
    F=np.eye(n,dtype=int)+2*V
    Sidx=[I.index(i) for i in S]
    rw=[int(np.abs(V[s]).sum()) for s in Sidx]
    print(dict(q=q,m=len(S),all_attained=ok,DF_eq_H=bool((D@F==Hm).all()),parity_bound=par_total,owner_total=sum(rw),
               rows=sorted(set(rw)),half_qm1=(q-1)/2,maxentry=int(np.abs(V).max())),flush=True)
