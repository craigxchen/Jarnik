import numpy as np
from common import *
rng=np.random.default_rng(1)
def F(v,t):
    n=v[:10]; x=v[10:15]; y=v[15:20]
    S=Svec(n); out=[x[i]**2+y[i]**2-S[i] for i in range(5)]
    for a,(i,j) in enumerate(E): out.append(x[i]*y[j]-x[j]*y[i]-t[a]*n[a])
    return np.array(out)
def jac(f,v,h=1e-7):
    f0=f(v); J=np.zeros((len(f0),len(v)))
    for k in range(len(v)):
        w=v.copy(); w[k]+=h; w2=v.copy(); w2[k]-=h; J[:,k]=(f(w)-f(w2))/(2*h)
    return J
for trial in range(6):
    if trial<3: t,n,P,th=real_sample(rng)
    else:
        th=rng.uniform(-np.pi,np.pi,5); th-=th[4]
        t,n,P,th=real_sample(rng,theta=th)
    v=np.concatenate([n,P.real,P.imag])
    r=F(v,t); 
    J=jac(lambda w:F(w,t),v)
    # scale columns
    sv=np.linalg.svd(J,compute_uv=False)
    rank=np.sum(sv>1e-8*sv[0])
    # tangent space: kernel of J
    U,s,Vt=np.linalg.svd(J); K=Vt[rank:].T  # 20 x (20-rank)
    rn=np.linalg.matrix_rank(K[:10],tol=1e-8*np.abs(K).max())
    # ratio angles: theta_i - theta_5 = arg(P_i)-arg(P_5)
    def ang(w):
        Pc=w[10:15]+1j*w[15:20]; a=np.angle(Pc); return (a[:4]-a[4])
    Ja=jac(ang,v); ra=np.linalg.matrix_rank(Ja@K,tol=1e-7)
    # their equations: c_i = 1/prod_j t_ij ; M c = 0
    M=np.array([[tval(t,i,j)*n[eid(i,j)] if i!=j else 0 for j in range(5)] for i in range(5)])
    c=np.array([1/np.prod([tval(t,i,j) for j in range(5) if j!=i]) for i in range(5)])
    # linear forms in n: row i coefficient on edge
    Lf=np.zeros((5,10))
    for i in range(5):
        for j in range(5):
            if j!=i: Lf[i,eid(i,j)]+=tval(t,i,j)*c[j]
    svM=np.linalg.svd(M,compute_uv=False)
    # sextic for triangle (0,1,2), l,m=3,4
    def A(i,j,k):
        l,m=[a for a in range(5) if a not in (i,j,k)]
        N=lambda a,b:n[eid(a,b)]; T=lambda a,b:t[eid(a,b)]
        return T(j,k)**2*N(j,k)*N(i,l)*N(i,m)+T(i,k)**2*N(i,k)*N(j,l)*N(j,m)-T(i,j)**2*N(i,j)*N(k,l)*N(k,m)
    worst=0
    for (i,j) in E:
        D=np.prod([n[eid(i,k)]*n[eid(j,k)] for k in range(5) if k not in (i,j)])
        for k in range(5):
            if k in (i,j): continue
            lhs=4*t[eid(i,k)]**2*t[eid(j,k)]**2*D-A(i,j,k)**2; rhs=4*(t[eid(i,j)]*t[eid(i,k)]*t[eid(j,k)])**2
            worst=max(worst,abs(lhs-rhs)/rhs)
    print(f"trial {trial}: t={t.astype(int)} resid={np.abs(r).max():.1e} rankJ={rank} dimX={20-rank} rank(dn)={rn} rank(dratio)={ra} "
          f"|Mc|={np.abs(M@c).max():.1e} rank(linear forms)={np.linalg.matrix_rank(Lf)} svM={svM[2]/svM[0]:.1e} sextic_relerr={worst:.1e}")
