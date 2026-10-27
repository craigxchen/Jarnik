# Exact minimum of sigma over F=0 ternary f with |supp f|>=4 (sdc_q), via HiGHS MILP.
# By the reduction theorem (hlight.md Thm 2.1) every other non-pair integral lambda has sigma >= (q-1)/2,
# so min(MILP value, (q-1)/2) is a lower bound for all non-pair lambda, and equals the global min when MILP < (q-1)/2
# (the global min over all non-pair lambda is then attained in this class).
# Symmetry: AGL(1,q) is sharply 2-transitive and preserves the patterns; negation preserves sigma:
#   WLOG f(0)=+1 and f(1)!=0.
import sys, time, numpy as np
from scipy.optimize import milp, LinearConstraint, Bounds
from hlib import chi_vec, sigma_formula
def solve(q,tl):
    ch=chi_vec(q)
    nx=q; N=3*q  # x(q) y(q) z(q)
    ix=lambda t:t; iy=lambda t:q+t; iz=lambda a:2*q+a
    rows=[];lo=[];hi=[]
    def add(d,l,h):
        r=np.zeros(N)
        for k,v in d.items(): r[k]+=v
        rows.append(r);lo.append(l);hi.append(h)
    for t in range(q): add({ix(t):1,iy(t):1},-np.inf,1)
    add({**{ix(t):1 for t in range(q)},**{iy(t):-1 for t in range(q)}},0,0)
    add({**{ix(t):1 for t in range(q)},**{iy(t):1 for t in range(q)}},4,np.inf)
    add({ix(0):1},1,1)
    add({ix(1):1,iy(1):1},1,1)
    for a in range(q):
        d={iz(a):1}
        for t in range(q):
            v=int(ch[(a-t)%q])
            if v: d[ix(t)]=d.get(ix(t),0)-v; d[iy(t)]=d.get(iy(t),0)+v
        add(dict(d),0,np.inf)
        d={iz(a):1}
        for t in range(q):
            v=int(ch[(a-t)%q])
            if v: d[ix(t)]=d.get(ix(t),0)+v; d[iy(t)]=d.get(iy(t),0)-v
        add(dict(d),0,np.inf)
    c=np.zeros(N); c[2*q:]=1
    lb=np.zeros(N); ub=np.ones(N); ub[2*q:]=np.inf
    integ=np.zeros(N); integ[:2*q]=1
    t0=time.time()
    res=milp(c,constraints=LinearConstraint(np.array(rows),lo,hi),bounds=Bounds(lb,ub),integrality=integ,
             options={"time_limit":tl,"disp":False,"mip_rel_gap":0})
    f=np.round(res.x[:q]-res.x[q:2*q]).astype(np.int64) if res.x is not None else None
    sig=sigma_formula(f,q,ch) if f is not None else None
    db=getattr(res,'mip_dual_bound',None)
    return res.status,sig,(None if db is None else db-q),f,time.time()-t0
if __name__=="__main__":
    tl=float(sys.argv[1])
    for q in [int(a) for a in sys.argv[2:]]:
        st,sig,dbs,f,tt=solve(q,tl)
        n=int(np.abs(f).sum()) if f is not None else None
        print(f"q={q} status={st} sigma*={sig} sigma*/r={sig/(2*q):.4f} dual_bound_sigma={dbs} (q-1)/2={(q-1)/2} n={n} supp={np.nonzero(f)[0].tolist()} vals={f[np.nonzero(f)[0]].tolist()} time={tt:.0f}s",flush=True)
