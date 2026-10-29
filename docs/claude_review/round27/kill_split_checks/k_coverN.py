# Rational curves of class N*x_F through e with 4N branches; F one-sided (each B in F hit N times, at distinct places).
# usage: k_coverN.py <signature> <N> <seconds> <seed>
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
roots=[B for B in F for _ in range(N)]          # 15N places
nr=len(roots); nb=4*N
A=[[j for j,B in enumerate(roots) if a in B and 4 not in B] for a in range(4)]
D=[[j for j,B in enumerate(roots) if 4 in B and a not in B] for a in range(4)]
fixed=np.array([0,1,-1],complex)
def unpack(z): return z[:nr], z[nr:nr+4], np.concatenate([fixed,z[nr+4:]])
def resid(z):
    r,c,S=unpack(z); out=[]
    for a in range(4):
        for si in S: out.append(c[a]*np.prod(si-r[A[a]])-np.prod(si-r[D[a]]))
    return np.array(out)
nz=nr+4+nb-3
def jac(z,h=1e-7):
    J=np.zeros((4*nb,nz),complex)
    for k in range(nz):
        dz=np.zeros(nz,complex); dz[k]=h; J[:,k]=(resid(z+dz)-resid(z-dz))/(2*h)
    return J
def newton(z,it=300):
    for t in range(it):
        f=resid(z); nf=np.linalg.norm(f)
        if nf<1e-12: return z,nf
        step=np.linalg.lstsq(jac(z),-f,rcond=None)[0]; lam=1.0
        while lam>1e-4:
            zn=z+lam*step
            if np.linalg.norm(resid(zn))<nf*(1-0.3*lam): break
            lam/=2
        z=zn
    return z,np.linalg.norm(resid(z))
def dominance(z,s=complex(0.37,0.21),h=1e-7):
    J=jac(z); _,sv,Vh=np.linalg.svd(J); K=Vh[J.shape[0]:].conj().T
    def Phi(zz,ss):
        r,c,_=unpack(zz); c=np.append(c,1.0)
        u=np.array([c[a]*np.prod([ss-r[j] for j,B in enumerate(roots) if a in B]) for a in range(5)])
        return np.log(u[:4]/u[4])
    cols=[(Phi(z+h*K[:,k],s)-Phi(z-h*K[:,k],s))/(2*h) for k in range(K.shape[1])]
    cols.append((Phi(z,s+h)-Phi(z,s-h))/(2*h))
    return np.linalg.svd(np.array(cols).T,compute_uv=False), sv[-1], K.shape[1]
print("profile",sorted(len(B) for B in F),"N=",N,"unknowns",nz,"equations",4*nb,flush=True)
t0=time.time(); tries=0; degen=0
while time.time()-t0<T:
    tries+=1
    z=(rng.normal(size=nz)+1j*rng.normal(size=nz))*10**rng.uniform(-0.5,0.7)
    z,nf=newton(z)
    if nf<1e-10:
        r,c,S=unpack(z)
        sep=min(min(abs(r[i]-r[j]) for i in range(nr) for j in range(i)), min(abs(si-rj) for si in S for rj in r), min(abs(S[i]-S[j]) for i in range(nb) for j in range(i)))
        if sep>1e-3 and min(abs(c))>1e-3 and min(abs(c-1))>1e-3:
            dsv,smin,kd=dominance(z)
            print(f"FOUND after {tries} tries: sep={sep:.3e}, kernel dim {kd}, eq-Jacobian sigma_min={smin:.3e}, dominance SVs={dsv}",flush=True)
            np.save(f"k_coverN_{N}_{int(sys.argv[4])}.npy",z); break
        degen+=1
print("tries",tries,"degenerate converged",degen,flush=True)
