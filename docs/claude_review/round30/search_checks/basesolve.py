# base-point form: roots q_k (any half plane), zero-sum sets F: sum_{k in F} Im q_k^m = 0, m=1..L
import numpy as np, itertools
from scipy.optimize import least_squares
def res_base(z,Fs,L):
    D=len(z)//2; q=z[:D]+1j*z[D:]
    P=np.array([np.imag(q**m) for m in range(1,L+1)])
    r=(Fs@P.T).ravel()
    return np.concatenate([r,[z[:D].sum(), (np.abs(q)**2).sum()-D]])
def degen(q):
    D=len(q); a=np.concatenate([q,np.conj(q)])
    d=min(abs(a[i]-a[j]) for i in range(2*D) for j in range(i+1,2*D))
    return d
def solve_base(Fs,L,starts=50,seed=0,y0=None):
    Fs=np.array(Fs,float); D=Fs.shape[1]; rng=np.random.default_rng(seed); out=[]
    for s in range(starts):
        z0=np.concatenate([rng.normal(0,1,D), (y0 if y0 is not None else rng.choice([-1,1],D))*np.abs(rng.normal(0.8,0.4,D))])
        r=least_squares(res_base,z0,args=(Fs,L),method=('lm' if Fs.shape[0]*L+2>=2*D else 'trf'),xtol=1e-15,ftol=1e-15,gtol=1e-15,max_nfev=5000)
        q=r.x[:D]+1j*r.x[D:]; cost=np.linalg.norm(r.fun)
        out.append((cost,degen(q),q))
    out.sort(key=lambda t:(t[0]>1e-10, -t[1]))
    return out
def zero_sum_sets(q,L,tol=1e-8):
    D=len(q); P=np.array([np.imag(q**m) for m in range(1,L+1)])
    sc=np.abs(P).max(axis=1,keepdims=True)
    S=[]
    for F in itertools.product([0,1],repeat=D):
        v=(P/sc)@np.array(F)
        if np.linalg.norm(v)<tol: S.append(F)
    return np.array(S)
