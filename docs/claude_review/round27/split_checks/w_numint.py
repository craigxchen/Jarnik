# Numerical intersection counts with W (evidence only): solve F(a,b,1,t) ∝ P(r) by Newton from random complex starts.
# mode 'line': P(r) = P0 + r P1 generic line -> #(W ∩ L) = deg W ; mode 'linee': line through e=(1,..,1) -> deg W - mult_e W
# mode 'curve': P(r) = lambda * C_R(r), a random rational curve with code profile x_R (R={2,3,4}) -> cost_{x_R}(g_W).
import numpy as np, sympy as sp, sys, time
mode=sys.argv[1]; nstart=int(sys.argv[2]); seed=int(sys.argv[3])
rng=np.random.default_rng(seed)
a,b,t,r=sp.symbols('a b t r')
c=1; d=a+b+c; m=a*b*c*d
S=a**2+b**2+c**2+d**2
za=(a**2*b+a**2*c+a*b**2+a*c**2-b**2*c-b*c**2)/2
den=t**2-m
xn=a*(t**2+m)-2*za*t          # x = xn/den
zn=za*den+t*(xn-a*den)        # z = zn/den
# homogenise F_j by den^4 (x, z have denominator den)
F0=(xn-a*den)*(xn-b*den)*(xn-c*den)*(xn+d*den)
F1=xn**4-(S/2)*xn**2*den**2+m*den**4+2*xn*zn*den**2
F2=(xn**2-b**2*den**2)*(xn**2-c**2*den**2)
F3=(xn**2-a**2*den**2)*(xn**2-c**2*den**2)
F4=(xn**2-a**2*den**2)*(xn**2-b**2*den**2)
F=[F0,F1,F2,F3,F4]
cplx=lambda n: rng.normal(size=n)+1j*rng.normal(size=n)
if mode in('line','linee'):
    P0=cplx(5) if mode=='line' else np.ones(5,dtype=complex)
    P1=cplx(5)
    P=[P0[j]+r*P1[j] for j in range(5)]
else:
    # code profile F_R, R={2,3,4} (0-based labels 0..4 = F_0..F_4): sets B with |B∩{2,3,4}| even
    import itertools
    FR=[set(B) for k in range(1,5) for B in itertools.combinations(range(5),k) if len(set(B)&{2,3,4})%2==0]
    assert len(FR)==15
    sB=cplx(15); lam=cplx(5)
    P=[lam[j]*sp.prod([(r-sB[i]) for i,B in enumerate(FR) if j in B]) for j in range(5)]
eqs=[sp.expand(F[j]*P[4]-F[4]*P[j]) for j in range(4)]
vars_=[a,b,t,r]
f=sp.lambdify(vars_,eqs,'numpy')
J=sp.lambdify(vars_,[[sp.diff(e,v) for v in vars_] for e in eqs],'numpy')
Fn=sp.lambdify(vars_,F+[den]+P,'numpy')
sols=[]
t0=time.time()
for it in range(nstart):
    v=cplx(4)*np.array([1,1,1,1.5])
    ok=False
    for k in range(80):
        try:
            fv=np.array(f(*v),dtype=complex); Jv=np.array(J(*v),dtype=complex)
            dv=np.linalg.solve(Jv,-fv)
        except Exception: break
        v=v+dv
        if np.linalg.norm(dv)<1e-12*(1+np.linalg.norm(v)): ok=True; break
        if np.linalg.norm(v)>1e8: break
    if not ok: continue
    vals=np.array(Fn(*v),dtype=complex)
    Fv,denv,Pv=vals[:5],vals[5],vals[6:]
    if abs(denv)<1e-6 or np.min(np.abs(Fv))<1e-8*np.max(np.abs(Fv)) or np.min(np.abs(Pv))<1e-8*np.max(np.abs(Pv)): continue
    rr=v[3]
    if not any(abs(rr-s)<1e-6*(1+abs(s)) for s in sols): sols.append(rr)
    if it%2000==0: print(it,len(sols),f"{time.time()-t0:.0f}s",flush=True)
print("mode",mode,"distinct intersection points found:",len(sols))
