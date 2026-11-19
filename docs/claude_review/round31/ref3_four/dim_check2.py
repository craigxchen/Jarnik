import numpy as np
from common import *
# exact-tangent version: tangent of X via analytic Jacobian, in log-scaled coords
rng=np.random.default_rng(1)
def Jan(n,P,t):
    x=P.real;y=P.imag;S=Svec(n)
    J=np.zeros((15,20))
    for i in range(5):
        J[i,10+i]=2*x[i]; J[i,15+i]=2*y[i]
        for j in range(5):
            if j!=i: J[i,eid(i,j)]=-S[i]/n[eid(i,j)]
    for a,(i,j) in enumerate(E):
        r=5+a; J[r,10+i]=y[j]; J[r,15+j]=x[i]; J[r,10+j]=-y[i]; J[r,15+i]=-x[j]; J[r,a]=-t[a]
    return J
for trial in range(12):
    th=rng.uniform(-np.pi,np.pi,5); th-=th[4]
    t,n,P,th=real_sample(rng,theta=th)
    J=Jan(n,P,t)
    # column scaling: n-columns by n, P columns by |P|
    sc=np.concatenate([n,np.abs(P),np.abs(P)]); Js=J*sc
    rs=Js/np.linalg.norm(Js,axis=1,keepdims=True)
    s=np.linalg.svd(rs,compute_uv=False)
    U,ss,Vt=np.linalg.svd(rs); K=(Vt[15:].T)  # scaled kernel
    sn=np.linalg.svd(K[:10],compute_uv=False)
    print(trial, "J sv tail", s[-2:].round(8)," dn sv",sn.round(8))
