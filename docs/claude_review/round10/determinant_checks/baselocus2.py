# EVIDENCE ONLY. Base locus Z_A = {delta : f_delta|_A in Poly_{<m}|_A} for small slab-boxes (full V(A,m), not V_H),
# and the tangent cone {u : <k,u>^m|_A in Poly_{<m}|_A}.  Grid search + local refinement, floating point.
import numpy as np, itertools, math, sys
from scipy.optimize import minimize
rng = np.random.default_rng(int(sys.argv[1]) if len(sys.argv)>1 else 1)
def slab_box(n,T,r,frac=0.5):
    d=n-1
    box=np.array(list(itertools.product(range(-T,T+1),repeat=d)))
    beta=rng.choice([-1,1],size=(r,n)).astype(float)
    X=box@(beta[:,:d]-beta[:,[d]]).T
    mx=np.abs(X).max(axis=1); w=np.sort(mx)[int(frac*len(mx))]
    return box[mx<=w].astype(float), beta, w
def polyQ(A,m):
    d=A.shape[1]
    mons=[a for t in range(m) for a in itertools.product(range(t+1),repeat=d) if sum(a)==t]
    P=np.stack([np.prod(A**np.array(a),axis=1) for a in mons],axis=1)
    U,s,_=np.linalg.svd(P,full_matrices=False); rk=int((s>s[0]*1e-10).sum())
    return U[:,:rk], rk
def res_batch(D,A,Q):
    F=np.exp(1j*(D@A.T)); R=F-(F@Q)@Q.T
    return (np.abs(R)**2).sum(axis=1)/A.shape[0]
def res_grad(x,A,Q):
    f=np.exp(1j*(A@x)); g=f-Q@(Q.T@f)
    return np.vdot(g,g).real/len(A), (2/len(A))*np.real((np.conj(g)*1j*f)@A)
def relation(x, vmax=6, tol=1e-6):
    d=len(x); full=np.append(x,0.0)
    best=None
    for v in itertools.product(range(-vmax,vmax+1),repeat=d):
        if all(t==0 for t in v): continue
        val=np.dot(v,x)/(2*math.pi)
        if abs(val-round(val))<tol:
            if best is None or sum(map(abs,v))<sum(map(abs,best)): best=v
    return best
for (n,T,r,m,gridN) in [(4,3,8,5,64),(4,4,8,6,64),(5,2,10,5,22)]:
    A,beta,w=slab_box(n,T,r); Q,rk=polyQ(A,m)
    print(f"\nn={n} T={T} r={r} m={m}: |A|={len(A)} rank Poly_<m|A={rk} dim V={len(A)-rk}",flush=True)
    d=n-1
    g=np.linspace(-math.pi,math.pi,gridN,endpoint=False)
    pts=np.array(list(itertools.product(g,repeat=d)))
    vals=np.concatenate([res_batch(pts[i:i+20000],A,Q) for i in range(0,len(pts),20000)])
    order=np.argsort(vals)[:400]
    zeros=[]
    for idx in order:
        x0=pts[idx]
        r_=minimize(res_grad,x0,args=(A,Q),jac=True,method='BFGS',options={'gtol':1e-13,'maxiter':400})
        if r_.fun<1e-18:
            x=(r_.x+math.pi)%(2*math.pi)-math.pi
            if not any(np.max(np.abs(((x-z)+math.pi)%(2*math.pi)-math.pi))<1e-5 for z in zeros): zeros.append(x)
    print(f"  grid min residual {vals.min():.2e}; distinct zeros after refinement: {len(zeros)}",flush=True)
    nondeg=0
    for z in zeros:
        rel=relation(z, vmax=4 if d==4 else 6)
        if rel is None: nondeg+=1
        print("   zero", np.round(z/(2*math.pi),6), "(units of 2pi); smallest integer relation:", rel, flush=True)
    print(f"  zeros with no integer relation |v|_inf<=6: {nondeg}",flush=True)
    # tangent cone: minimise tau(u) = ||(I-P)<k,u>^m|| / ||<k,u>^m|| over unit u
    def tau(u):
        u=u/np.linalg.norm(u); f=(A@u)**m; gg=f-Q@(Q.T@f); return np.dot(gg,gg)/np.dot(f,f)
    best=[]
    for s in range(300):
        rr=minimize(tau,rng.normal(size=d),method='Nelder-Mead',options={'xatol':1e-12,'fatol':1e-20,'maxiter':3000})
        best.append((rr.fun, rr.x/np.linalg.norm(rr.x)))
    best.sort(key=lambda t:t[0])
    print(f"  tangent cone: min tau = {best[0][0]:.3e} at u = {np.round(best[0][1],4)}; 5th smallest {best[4][0]:.3e}",flush=True)
