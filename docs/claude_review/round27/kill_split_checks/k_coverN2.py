# Faster N-fold search (analytic Jacobian, ratio-form residuals g_{a,i} - 1, g = c_a prod_A(s_i-r)/prod_D(s_i-r)).
# usage: k_coverN2.py <signature> <N> <seconds> <seed>
import numpy as np, itertools, sys, time
full=frozenset(range(5))
Bs=[frozenset(s) for r in range(1,5) for s in itertools.combinations(range(5),r)]
cuts=[B for B in Bs if 0 in B]
sig=eval(sys.argv[1]); N=int(sys.argv[2]); T=float(sys.argv[3]); rng=np.random.default_rng(int(sys.argv[4]))
F=None
for choice in range(1<<15):
    FF=[cuts[i] if (choice>>i)&1 else full-cuts[i] for i in range(15)]
    X=[sum(1 for B in FF if a in B) for a in range(5)]
    if len(set(X))==1 and sorted(len(B) for B in FF)==sig: F=FF; break
roots=[B for B in F for _ in range(N)]; nr=len(roots); nb=4*N; nfree=nb-3
SA=np.zeros((4,nr)); SD=np.zeros((4,nr))
for a in range(4):
    for j,B in enumerate(roots):
        if a in B and 4 not in B: SA[a,j]=1
        if 4 in B and a not in B: SD[a,j]=1
W=SA-SD   # exponent of (s_i - r_j) in g_a
fixed=np.array([0,1,-1],complex)
nz=nr+4+nfree
def unpack(z): return z[:nr], z[nr:nr+4], np.concatenate([fixed,z[nr+4:]])
Aidx=[np.where(SA[a]>0)[0] for a in range(4)]; Didx=[np.where(SD[a]>0)[0] for a in range(4)]
def prods_excl(v):
    # products of all entries except one, without division
    n=len(v); pre=np.ones(n+1,complex); suf=np.ones(n+1,complex)
    for k in range(n): pre[k+1]=pre[k]*v[k]
    for k in range(n-1,-1,-1): suf[k]=suf[k+1]*v[k]
    return pre[n], np.array([pre[k]*suf[k+1] for k in range(n)])
def F_J(z):
    r,c,S=unpack(z)
    res=np.zeros(4*nb,complex); J=np.zeros((4*nb,nz),complex)
    for a in range(4):
        for i in range(nb):
            row=a*nb+i
            dA=S[i]-r[Aidx[a]]; dD=S[i]-r[Didx[a]]
            pA,eA=prods_excl(dA); pD,eD=prods_excl(dD)
            res[row]=c[a]*pA-pD
            J[row,Aidx[a]]+=-c[a]*eA
            J[row,Didx[a]]+=eD
            J[row,nr+a]=pA
            if i>=3: J[row,nr+4+i-3]=c[a]*eA.sum()-eD.sum()
    return res,J
def newton(z,it=200):
    for t in range(it):
        f,J=F_J(z); nf=np.linalg.norm(f)
        if not np.isfinite(nf): return z,np.inf
        if nf<1e-12: return z,nf
        step=np.linalg.lstsq(J,-f,rcond=None)[0]; lam=1.0
        while lam>1e-4:
            zn=z+lam*step; fn,_=F_J(zn)
            if np.isfinite(np.linalg.norm(fn)) and np.linalg.norm(fn)<nf*(1-0.3*lam): break
            lam/=2
        z=zn
    f,_=F_J(z); return z,np.linalg.norm(f)
def Phi(z,s):
    r,c,_=unpack(z); c=np.append(c,1.0)
    u=np.array([c[a]*np.prod([s-r[j] for j,B in enumerate(roots) if a in B]) for a in range(5)])
    return np.log(u[:4]/u[4])
def dominance(z,s=complex(0.37,0.21),h=1e-7):
    _,J=F_J(z); U,sv,Vh=np.linalg.svd(J); K=Vh[J.shape[0]:].conj().T
    cols=[(Phi(z+h*K[:,k],s)-Phi(z-h*K[:,k],s))/(2*h) for k in range(K.shape[1])]
    cols.append((Phi(z,s+h)-Phi(z,s-h))/(2*h))
    return np.linalg.svd(np.array(cols).T,compute_uv=False), sv[-1], K.shape[1]
print("profile",sorted(len(B) for B in F),"N=",N,"unknowns",nz,"equations",4*nb,flush=True)
t0=time.time(); tries=0; degen=0; best=np.inf
while time.time()-t0<T:
    tries+=1
    z=(rng.normal(size=nz)+1j*rng.normal(size=nz))*1.5
    try:
        z,nf=newton(z)
    except Exception:
        continue
    best=min(best,nf)
    if nf<1e-10:
        r,c,S=unpack(z)
        sep=min(min(abs(r[i]-r[j]) for i in range(nr) for j in range(i)), min(abs(si-rj) for si in S for rj in r), min(abs(S[i]-S[j]) for i in range(nb) for j in range(i)))
        if sep>1e-3 and min(abs(c))>1e-3 and min(abs(c-1))>1e-3 and max(abs(r))<1e3 and max(abs(c))<1e6 and max(abs(S))<1e3:
            dsv,smin,kd=dominance(z)
            if min(dsv)<1e-6 or smin<1e-8: degen+=1; continue
            print(f"FOUND after {tries} tries ({time.time()-t0:.0f}s): sep={sep:.3e}, kernel dim {kd}, eq-Jacobian sigma_min={smin:.3e}, dominance SVs={dsv}",flush=True)
            np.save(f"k_coverN2_{N}_{int(sys.argv[4])}.npy",z); break
        degen+=1
print("tries",tries,"degenerate converged",degen,"best residual",best,flush=True)
