# Global min of sigma_lambda over integral sum-zero lambda with support >= 3 (non-pair), sdc_q, via HiGHS MILP.
# sigma = sum_S |lambda(S)| - q. Translation symmetry: force lambda_0 != 0 (support meets F_q) -- WLOG.
import sys, numpy as np, time
from scipy.optimize import milp, LinearConstraint, Bounds
sys.path.insert(0,'/tmp/claude-0/-home-user-Jarnik/ba4929a4-8665-5790-985a-3f09130a14b9/scratchpad/round16/referee_barrier_checks')
from sdc_sym import check
def solve(q,B=3,tl=600,minsupp=3):
    _,_,_,M=check(q); n=q+1; ns=2*q
    # vars: p(n), m(n), s(n) bin sign, b(n) bin support, t(ns)
    N=4*n+ns
    ip=lambda x:x; im=lambda x:n+x; isg=lambda x:2*n+x; ib=lambda x:3*n+x; it=lambda S:4*n+S
    c=np.zeros(N); c[4*n:]=1
    A=[];lo=[];hi=[]
    def row(d,l,h):
        r=np.zeros(N)
        for k,v in d.items(): r[k]+=v
        A.append(r);lo.append(l);hi.append(h)
    row({**{ip(x):1 for x in range(n)},**{im(x):-1 for x in range(n)}},0,0)
    for x in range(n):
        row({ip(x):1,isg(x):-B},-np.inf,0)
        row({im(x):1,isg(x):B},-np.inf,B)
        row({ip(x):1,im(x):1,ib(x):-1},0,np.inf)
        row({ip(x):1,im(x):1,ib(x):-B},-np.inf,0)
    row({ib(x):1 for x in range(n)},minsupp,np.inf)
    row({ib(0):1},1,1)
    for S in range(ns):
        d={it(S):1}
        for x in range(n):
            if M[S,x]: d[ip(x)]=d.get(ip(x),0)-1; d[im(x)]=d.get(im(x),0)+1
        row(dict(d),0,np.inf)
        d={it(S):1}
        for x in range(n):
            if M[S,x]: d[ip(x)]=d.get(ip(x),0)+1; d[im(x)]=d.get(im(x),0)-1
        row(dict(d),0,np.inf)
    lb=np.zeros(N); ub=np.full(N,np.inf)
    ub[:2*n]=B; ub[2*n:4*n]=1
    integ=np.ones(N); integ[4*n:]=0
    t0=time.time()
    res=milp(c,constraints=LinearConstraint(np.array(A),lo,hi),bounds=Bounds(lb,ub),integrality=integ,options={"time_limit":tl,"disp":False})
    lam=np.round(res.x[:n]-res.x[n:2*n]).astype(int)
    sig=int(np.abs(M@lam).sum())-q
    return res.status,res.fun-q if res.fun is not None else None,sig,lam,getattr(res,'mip_dual_bound',None),time.time()-t0
if __name__=="__main__":
    for q in [int(a) for a in sys.argv[1:]]:
        st,f,sig,lam,db,tt=solve(q)
        print(f"q={q} status={st} sigma={sig} sigma/r={sig/(2*q):.4f} dual_bound_sigma={None if db is None else db-q} norm1={np.abs(lam).sum()} supp={np.count_nonzero(lam)} lam={lam.tolist()} time={tt:.0f}s",flush=True)
