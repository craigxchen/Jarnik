import json, sys, numpy as np, itertools, mpmath as mp
from scipy.optimize import least_squares
from collections import Counter
exec(open('equi_num.py').read().split('for ci,cl in enumerate(cls):')[0].replace("D=int(sys.argv[1]); M=int(sys.argv[2]); starts=int(sys.argv[3]) if len(sys.argv)>3 else 200","D=8;M=5;starts=0"))
mp.mp.dps=50
def polish(x,F,its=30):
    # Newton (Gauss-Newton) in mpmath on the square-ish system: use least squares steps with numpy then mp refine via pinv
    x=[mp.mpf(v) for v in x]
    for it in range(its):
        rho=[mp.mpc(x[k],x[D+k]) for k in range(D)]
        def R(xx):
            rr=[mp.mpc(xx[k],xx[D+k]) for k in range(D)]
            out=[]
            for i in range(F.shape[0]):
                for j in range(1,L+1):
                    out.append(sum(F[i,k]*mp.im(rr[k]**j) for k in range(D)))
            out.append(sum(xx[:D])); out.append(sum(abs(r)**2 for r in rr)-D)
            return out
        r0=R(x)
        h=mp.mpf('1e-30')
        J=mp.matrix(len(r0),2*D)
        for v in range(2*D):
            xp=list(x); xp[v]+=h
            rp=R(xp)
            for u in range(len(r0)): J[u,v]=(rp[u]-r0[u])/h
        # minimal-norm step
        Jt=J.T
        step=Jt*mp.inverse(J*Jt)*mp.matrix(r0)
        x=[x[v]-step[v] for v in range(2*D)]
    return x, max(abs(v) for v in R(x))
rng=np.random.default_rng(5)
for ci in (0,1):
    F=build(cls[ci])
    for s in range(300):
        sol=least_squares(resid,rng.normal(size=2*D),args=(F,),xtol=1e-15,ftol=1e-15,gtol=1e-15,max_nfev=4000)
        if np.linalg.norm(sol.fun)<1e-10 and nondeg(sol.x)>0.03: break
    x,res=polish(sol.x,F)
    rho=[mp.mpc(x[k],x[D+k]) for k in range(D)]
    print('class',ci,'types',cls[ci]); print(' residual',mp.nstr(res,5),'nondeg',nondeg(np.array([float(v) for v in x])))
    print(' roots',[mp.nstr(r,8) for r in rho])
    # full fibre count over all 256 points
    C=Counter()
    for sg in itertools.product([1,-1],repeat=D):
        key=tuple(mp.nstr(sum(sg[k]*mp.im(rho[k]**j) for k in range(D)),20) for j in range(1,L+1))
        C[key]+=1
    print(' fibre sizes top',sorted(C.values())[-4:])
