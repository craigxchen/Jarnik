# Numerical intersection counts with W (EVIDENCE ONLY; random-start Newton can miss solutions, so counts are lower bounds).
# Unknowns (x,a,b,z,r), c = 1. W-point Phi(a,b,1,x,z) with z^2 = m x^2 + B. Target curve P(r) in P^4.
# Equations: z^2 - m x^2 - B = 0, and F_j(x,a,b,z)*P_4(r) - F_4*P_j(r) = 0 (j = 0..3).
# modes: 'line'  P = P0 + r P1 (generic)            -> deg W
#        'linee' P = (1,..,1) + r v (through e)      -> deg W - mult_e W
#        'curve' P = lam * C_R(r), code profile x_R  -> cost_{x_R}(g_W)
import numpy as np, sys, itertools, time
mode=sys.argv[1]; NST=int(sys.argv[2]); seed=int(sys.argv[3])
rng=np.random.default_rng(seed)
cplx=lambda n: rng.normal(size=n)+1j*rng.normal(size=n)
if mode=='line':
    P0,P1=cplx(5),cplx(5); Pf=lambda r: P0+r*P1
elif mode=='linee':
    v=cplx(5); Pf=lambda r: np.ones(5)+r*v
else:
    FR=[set(B) for k in range(1,5) for B in itertools.combinations(range(5),k) if len(set(B)&{2,3,4})%2==0]
    sB=cplx(15); lam=cplx(5)
    def Pf(r):
        return np.array([lam[j]*np.prod([(r-sB[i]) for i,B in enumerate(FR) if j in B]) for j in range(5)])
def eqs(w):
    x,a,b,z,r=w; c=1.0; d=a+b+c; m=a*b*c*d; sq=[a*a,b*b,c*c,d*d]; S=sum(sq)
    e3=sq[0]*sq[1]*sq[2]+sq[0]*sq[1]*sq[3]+sq[0]*sq[2]*sq[3]+sq[1]*sq[2]*sq[3]; B=(e3-S*m)/4
    F=np.array([(x-a)*(x-b)*(x-c)*(x+d), x**4-S/2*x*x+m+2*x*z, (x*x-b*b)*(x*x-c*c),(x*x-a*a)*(x*x-c*c),(x*x-a*a)*(x*x-b*b)])
    P=Pf(r)
    return np.concatenate([[z*z-m*x*x-B],F[:4]*P[4]-F[4]*P[:4]]), F, P
def J(w,h=1e-7):
    f0=eqs(w)[0]; M=np.zeros((5,5),dtype=complex)
    for i in range(5):
        e=np.zeros(5,dtype=complex); e[i]=h*(1+abs(w[i])); M[:,i]=(eqs(w+e)[0]-f0)/e[i]
    return M
sols=[]; t0=time.time(); hist=[]
with np.errstate(all='ignore'):
    for it in range(NST):
        w=cplx(5)*np.array([2,2,2,4,2])
        ok=False
        for k in range(60):
            f=eqs(w)[0]
            try: dw=np.linalg.solve(J(w),-f)
            except Exception: break
            w=w+dw
            if not np.all(np.isfinite(w)) or np.abs(w).max()>1e7: break
            if np.linalg.norm(dw)<1e-12*(1+np.linalg.norm(w)): ok=True; break
        if not ok: continue
        f,F,P=eqs(w)
        if np.abs(f).max()>1e-8*(1+np.abs(F).max()*np.abs(P).max()): continue
        if np.abs(F).min()<1e-7*np.abs(F).max() or np.abs(P).min()<1e-7*np.abs(P).max(): continue
        if mode=='linee' and abs(w[4])<1e-6: continue
        r=w[4]
        if not any(abs(r-s)<1e-7*(1+abs(s)) for s in sols): sols.append(r)
        if (it+1)%(NST//10)==0: hist.append(len(sols)); print(it+1,len(sols),f"{time.time()-t0:.0f}s",flush=True)
print("mode",mode,"seed",seed,"distinct intersection points found:",len(sols),"history",hist)
