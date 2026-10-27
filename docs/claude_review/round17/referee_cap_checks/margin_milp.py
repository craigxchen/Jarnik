# Referee: HiGHS MILP check (evidence for bounded entries) of the input of Prop 6.1:
# sum_P |v(P)| >= q+1 for every nonzero integral v (|v_x| <= B), over the 2q sdc_q relative patterns.
# WLOG v_0 >= 1 (translation invariance of the pattern set, and v -> -v).
import numpy as np, scipy.optimize as so, sys
def chi(a,q):
    a%=q
    return 0 if a==0 else (1 if pow(a,(q-1)//2,q)==1 else -1)
for q in [7,11,19,23,31,43]:
    for B in [3]:
        pats=[[a]+[x for x in range(q) if chi(a-x,q)==1] for a in range(q)]+[[a]+[x for x in range(q) if chi(a-x,q)==-1] for a in range(q)]
        np_=2*q; nv=q
        # vars: v (q, integer in [-B,B]), t (2q, continuous >=0)
        c=np.concatenate([np.zeros(nv),np.ones(np_)])
        rows=[];lb=[];ub=[]
        for i,S in enumerate(pats):
            r=np.zeros(nv+np_); r[S]=1; r[nv+i]=-1; rows.append(r); lb.append(-np.inf); ub.append(0)   # v(P)-t<=0
            r=np.zeros(nv+np_); r[S]=-1; r[nv+i]=-1; rows.append(r); lb.append(-np.inf); ub.append(0) # -v(P)-t<=0
        cons=so.LinearConstraint(np.array(rows),lb,ub)
        lo=np.concatenate([-B*np.ones(nv),np.zeros(np_)]); hi=np.concatenate([B*np.ones(nv),np.full(np_,np.inf)])
        lo[0]=1
        integ=np.concatenate([np.ones(nv),np.zeros(np_)])
        res=so.milp(c,constraints=cons,bounds=so.Bounds(lo,hi),integrality=integ,options={"time_limit":600})
        v=np.round(res.x[:nv]).astype(int)
        print(f"q={q} B={B}: status={res.status} ({res.message}) min sum_P|v(P)| = {res.fun:.6g}  (q+1={q+1})  argmin v={v.tolist()}",flush=True)
