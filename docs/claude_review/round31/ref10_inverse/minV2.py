import sys, math, random
import numpy as np
from minV import setup, solve_col, run, chi
from scipy.optimize import milp, LinearConstraint, Bounds
from sympy import nextprime
# (a) random S of density rho in the class, weights log p_k with split primes p_k in [P,2P]
def weights(n,P,seed):
    random.seed(seed); ps=set()
    while len(ps)<n:
        p=nextprime(random.randint(P,2*P))
        if p%4==1 and p<=2*P: ps.add(p)
    return np.array([math.log(p) for p in sorted(ps)])
# (b) free targets: min sum_{i in S}|v_i| + sum_rows |e_r| s.t. D v - e = b
def free_target(q):
    I,S,i0,rows,D,Hm,Bm=setup(q); n=len(I); r=len(rows)
    tot=0
    for t in range(n):
        b=Bm[:,t]
        # vars v(n) tv(n) e(r) te(r)
        N=2*n+2*r
        c=np.zeros(N)
        for j,k in enumerate(I): c[n+j]=1.0 if k in S else 1e-4
        c[2*n+r:]=1.0
        A=np.zeros((r,N)); A[:,:n]=D; A[:,2*n:2*n+r]=-np.eye(r)
        cons=[LinearConstraint(A,b,b)]
        for (o,w) in ((0,n),(2*n,r)):
            M1=np.zeros((w,N)); M1[:,o:o+w]=np.eye(w); M1[:,o+w:o+2*w]=-np.eye(w)
            M2=np.zeros((w,N)); M2[:,o:o+w]=-np.eye(w); M2[:,o+w:o+2*w]=-np.eye(w)
            cons+= [LinearConstraint(M1,-np.inf,0),LinearConstraint(M2,-np.inf,0)]
        integ=np.zeros(N); integ[:n]=1; integ[2*n:2*n+r]=1
        lb=np.full(N,-60.0); ub=np.full(N,60.0)
        res=milp(c,constraints=cons,integrality=integ,bounds=Bounds(lb,ub),options={'time_limit':120})
        assert res.status==0
        tot+=round(res.fun)
    return tot
if __name__=="__main__":
    mode=sys.argv[1]
    if mode=="rand":
        for q in [19,23,31]:
            I,S,*_=setup(q)
            for seed in range(2):
                random.seed(100+seed); Ssub=sorted(random.sample(S,max(3,int(0.6*len(S)))))
                w=weights(q-1,10**4,seed)
                R=run(q,S=Ssub,wts=w); R.pop('rows')
                print(R,flush=True)
    elif mode=="all":
        for q in [11,19,23]:
            R=run(q,owner_only=False); R.pop('rows'); print(R,flush=True)
    elif mode=="free":
        for q in [11,19,23]:
            I,S,*_=setup(q); print(q,'free-target min owner-twist+target excess =',free_target(q),' (q-3)/2*m =',(q-3)//2*len(S),flush=True)
