# Exact column minimum of sum_{owners}|v_i| over ALL integer V (non-owners eliminated via Smith form, unbounded)
import sys, numpy as np
from sympy import Matrix, ZZ
from sympy.matrices.normalforms import smith_normal_decomp
from minV import setup
from scipy.optimize import milp, LinearConstraint, Bounds
for q in map(int,sys.argv[1:]):
    I,S,i0,rows,D,Hm,Bm=setup(q); n=len(I)
    own=[t for t,k in enumerate(I) if k in S]; non=[t for t,k in enumerate(I) if k not in S]
    DN=Matrix(D[:,non].tolist()); DS=D[:,own].astype(float)
    Sm,Pm,Qm=smith_normal_decomp(DN,domain=ZZ); assert Sm==Pm*DN*Qm
    r=min(Sm.shape); dd=[int(Sm[i,i]) for i in range(r)]; rank=sum(1 for x in dd if x)
    Pn=[[int(x) for x in Pm.row(i)] for i in range(Pm.rows)]
    conds=[(Pn[i],abs(dd[i])) for i in range(rank) if abs(dd[i])!=1]+[(Pn[i],0) for i in range(rank,Sm.rows)]
    Bf={i0:np.zeros(n,dtype=int)}
    for r_,i in enumerate(rows): Bf[i]=Bm[r_]
    no=len(own); nc=len(conds); BND=7
    tot=0; partot=0; notopt=0; excess=[]; rowsum=np.zeros(no,dtype=int)
    for t in range(n):
        b=Bm[:,t]
        # vars x(no) a(no) z(nc)
        N=2*no+nc; c=np.concatenate([np.zeros(no),np.ones(no),np.zeros(nc)])
        cons=[]
        if nc:
            A=np.zeros((nc,N)); rhs=np.zeros(nc)
            for j,(w,d) in enumerate(conds):
                w=np.array(w,dtype=float); A[j,:no]=w@DS; A[j,2*no+j]=d; rhs[j]=w@b
            cons.append(LinearConstraint(A,rhs,rhs))
        E=np.eye(no); Z=np.zeros((no,nc))
        cons+=[LinearConstraint(np.hstack([E,-E,Z]),-np.inf,0),LinearConstraint(np.hstack([-E,-E,Z]),-np.inf,0)]
        lb=np.concatenate([-BND*np.ones(no),np.zeros(no),-np.inf*np.ones(nc)]); ub=np.concatenate([BND*np.ones(no),BND*np.ones(no),np.inf*np.ones(nc)])
        res=milp(c,constraints=cons,integrality=np.concatenate([np.ones(no),np.zeros(no),np.ones(nc)]),bounds=Bounds(lb,ub),options={'time_limit':60})
        pm=min(sum((cc+Bf[I[o]][t])%2 for o in own) for cc in (0,1)); partot+=pm
        if res.x is None: print('col',t,'no solution found',res.message,flush=True); notopt+=1; continue
        x=np.round(res.x[:no]).astype(int)
        # exact integrality check of residual solvability: verify conditions with integers
        rhsr=b-D[:,own]@x
        okc=all((sum(int(wi)*int(ri) for wi,ri in zip(w,rhsr)) % d==0) if d else sum(int(wi)*int(ri) for wi,ri in zip(w,rhsr))==0 for w,d in conds)
        assert okc
        val=int(np.abs(x).sum()); tot+=val; rowsum+=np.abs(x)
        if res.status!=0: notopt+=1
        print(' col',t,'parity',pm,'val',val,'status',res.status,flush=True)
        if val>pm: excess.append((t,pm,val,res.status))
    print(dict(q=q,m=no,invariants_DN=sorted(set(abs(x) for x in dd)),parity_bound=partot,ILP_total=tot,
               cols_not_proven_opt=notopt,cols_above_parity=len(excess),excess=excess[:6],
               owner_rows=sorted(set(rowsum.tolist())),row_target=(q-3)/2),flush=True)
