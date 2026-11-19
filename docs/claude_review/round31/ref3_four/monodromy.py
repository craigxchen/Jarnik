# Degree of V_t by monodromy, using ONLY the star-vector incidence (no use of Astra's equations).
# unknowns z=(x1..x5,y1..y4), y5=0; n_ij=(x_i y_j-x_j y_i)/t_ij ; eqs x_i^2+y_i^2=prod_j n_ij (5), Lam n = c (4)
import numpy as np, sys
from common import *
seed=int(sys.argv[1]) if len(sys.argv)>1 else 0
rng=np.random.default_rng(seed)
th=rng.uniform(-np.pi,np.pi,5); th-=th[4]
t,n0,P0,_=real_sample(rng,theta=th)
scale=np.abs(P0).mean(); 
Lam=rng.normal(size=(4,10))+1j*rng.normal(size=(4,10))
def nvec(z):
    x=z[:5]; y=np.concatenate([z[5:9],[0]])
    return np.array([(x[i]*y[j]-x[j]*y[i])/t[a] for a,(i,j) in enumerate(E)])
def Fz(z,c):
    x=z[:5]; y=np.concatenate([z[5:9],[0]]); n=nvec(z); S=Svec(n)
    return np.concatenate([x**2+y**2-S, Lam@n-c])
def Jz(z):
    x=z[:5]; y=np.concatenate([z[5:9],[0]]); n=nvec(z); S=Svec(n)
    dn=np.zeros((10,9),complex)
    for a,(i,j) in enumerate(E):
        dn[a,i]+= y[j]/t[a]; dn[a,j]+= -y[i]/t[a]
        if j<4: dn[a,5+j]+= x[i]/t[a]
        if i<4: dn[a,5+i]+= -x[j]/t[a]
    J=np.zeros((9,9),complex)
    for i in range(5):
        J[i,i]+=2*x[i]
        if i<4: J[i,5+i]+=2*y[i]
        for j in range(5):
            if j!=i: J[i]-= np.prod([n[eid(i,k)] for k in range(5) if k not in (i,j)])*dn[eid(i,j)]
    J[5:]=Lam@dn
    return J
z0=np.concatenate([P0.real,P0.imag[:4]]).astype(complex)
c0=Lam@nvec(z0)
assert np.abs(Fz(z0,c0)).max()<1e-8*np.abs(c0).max()+1e-6
def track(z,ca,cb,gam):
    try:
        return track0(z,ca,cb)
    except Exception:
        return None
def track0(z,ca,cb):
    # path c(s)=(1-s)ca+s cb along straight segment; s in [0,1]
    s=0.0; h=0.02
    while s<1:
        h=min(h,1-s)
        dc=cb-ca
        # predictor RK4 on dz/ds = -J^{-1} dc  (dF/dc=-I on last 4)
        def vel(zz):
            rhs=np.concatenate([np.zeros(5),dc]); return np.linalg.solve(Jz(zz),rhs)
        k1=vel(z); k2=vel(z+h/2*k1); k3=vel(z+h/2*k2); k4=vel(z+h*k3)
        zp=z+h/6*(k1+2*k2+2*k3+k4)
        cp=ca+(s+h)*dc
        ok=False
        for it in range(6):
            dz=np.linalg.solve(Jz(zp),-Fz(zp,cp)); zp=zp+dz
            if np.linalg.norm(dz)<1e-10*(1+np.linalg.norm(zp)): ok=True; break
        if ok and np.linalg.norm(zp-z)<0.3*(1+np.linalg.norm(z)):
            z=zp; s+=h; h=min(h*1.6,0.1)
        else:
            h/=2
            if h<1e-9 or not np.all(np.isfinite(z)) or np.linalg.norm(z)>1e8*scale: return None
    return z
sols=[z0]
def key(z): return z
def isnew(z):
    return all(np.linalg.norm(z-w)>1e-6*(1+np.linalg.norm(w)) for w in sols)
stall=0; loops=0
import time; T0=time.time()
while stall<6 and loops<80:
    loops+=1
    c1=c0+np.abs(c0).max()*(rng.normal(size=4)+1j*rng.normal(size=4))
    c2=c0+np.abs(c0).max()*(rng.normal(size=4)+1j*rng.normal(size=4))
    added=False
    for z in list(sols):
        w=track(z,c0,c1,None)
        if w is None: continue
        w=track(w,c1,c2,None)
        if w is None: continue
        w=track(w,c2,c0,None)
        if w is None: continue
        if np.abs(Fz(w,c0)).max()<1e-6*np.abs(c0).max() and isnew(w):
            sols.append(w); added=True
    stall=0 if added else stall+1
    print('loop',loops,'sols',len(sols),'time',round(time.time()-T0),flush=True)
ns=[nvec(z) for z in sols]
distinct=[]
for v in ns:
    if all(np.linalg.norm(v-u)>1e-6*np.linalg.norm(u) for u in distinct): distinct.append(v)
print("seed",seed,"t",t.astype(int),"loops",loops,"#star solutions",len(sols),"#distinct n",len(distinct))
np.save(f"mono_n_{seed}.npy",np.array(distinct)); np.save(f"mono_t_{seed}.npy",t)
