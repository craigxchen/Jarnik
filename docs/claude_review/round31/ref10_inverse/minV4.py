# Exact: non-owner unknowns eliminated by Smith form of D_N (no bounds); owners in {-1,0,1} with optimal parity
import sys, numpy as np
from sympy import Matrix, ZZ
from sympy.matrices.normalforms import smith_normal_decomp
from minV import setup
from scipy.optimize import milp, LinearConstraint, Bounds
def smith(M):
    Sm,Pm,Qm=smith_normal_decomp(M,domain=ZZ)  # Sm = Pm*M*Qm
    return Sm,Pm,Qm
for q in map(int,sys.argv[1:]):
    I,S,i0,rows,D,Hm,Bm=setup(q); n=len(I)
    own=[t for t,k in enumerate(I) if k in S]; non=[t for t,k in enumerate(I) if k not in S]
    DN=Matrix(D[:,non].tolist()); DS=D[:,own]
    Sm,Pm,Qm=smith(DN)
    assert Sm==Pm*DN*Qm
    r=min(Sm.shape); dd=[Sm[i,i] for i in range(r)]
    rank=sum(1 for x in dd if x!=0)
    print('q',q,'m',len(S),'D_N shape',DN.shape,'invariants',sorted(set(abs(x) for x in dd)),'rank',rank,flush=True)
    Pn=np.array(Pm.tolist(),dtype=object)
    # conditions: (Pn @ rhs)_i == 0 mod d_i for i<rank (d_i!=1), and ==0 for i>=rank
    conds=[(Pn[i],int(dd[i])) for i in range(rank) if abs(dd[i])!=1]+[(Pn[i],0) for i in range(rank,Sm.shape[0])]
    Bf={i0:np.zeros(n,dtype=int)}
    for r_,i in enumerate(rows): Bf[i]=Bm[r_]
    tot=0; partot=0; worst=[]
    for t in range(n):
        b=Bm[:,t]
        best=None
        for c in (0,1):
            par=[(c+Bf[I[o]][t])%2 for o in own]; cnt=sum(par)
            # vars: x (owner v, |own|), z (len conds); constraint w.(b - DS x) = d z
            no=len(own); nc=len(conds)
            if nc==0:
                sol=cnt
            else:
                A=np.zeros((nc,no+nc)); rhs=np.zeros(nc)
                for j,(w,d) in enumerate(conds):
                    w=np.array([int(x) for x in w],dtype=float)
                    A[j,:no]=-(w@DS); A[j,no+j]=-d; rhs[j]=-(w@b)
                # parity: owner even -> 0, odd -> +-1 : x = par + 2y? simpler: x in {-1,0,1}, x even ->0
                lb=np.array([0 if p==0 else -1 for p in par]+[-1e9]*nc,dtype=float)
                ub=np.array([0 if p==0 else 1 for p in par]+[1e9]*nc,dtype=float)
                # forbid x=0 on odd: x = 2s-1 with s binary
                Am=np.zeros((no,no+nc+no)); 
                A2=np.hstack([A,np.zeros((nc,no))])
                Ap=np.zeros((no,no+nc+no)); lo=np.zeros(no); hi=np.zeros(no)
                for o in range(no):
                    Ap[o,o]=1
                    if par[o]==1: Ap[o,no+nc+o]=-2; lo[o]=hi[o]=-1
                lbs=np.concatenate([lb,np.zeros(no)]); ubs=np.concatenate([ub,np.ones(no)])
                res=milp(np.zeros(no+nc+no),constraints=[LinearConstraint(A2,rhs,rhs),LinearConstraint(Ap,lo,hi)],
                         integrality=np.ones(no+nc+no),bounds=Bounds(lbs,ubs),options={'time_limit':120})
                sol=cnt if res.x is not None else None
            if sol is not None and (best is None or sol<best): best=sol
        par_min=min(sum((c+Bf[I[o]][t])%2 for o in own) for c in (0,1))
        partot+=par_min
        if best is None or best>par_min: worst.append((t,par_min,best))
        tot+= best if best is not None else 10**6
    print(' parity bound',partot,' attained with owner entries in {-1,0,1}:',tot,' columns not attained:',worst,flush=True)
