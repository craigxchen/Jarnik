import numpy as np, sys
from common import *
s=sys.argv[1]
N=np.load(f"mono_n_{s}.npy"); t=np.load(f"mono_t_{s}.npy")
c=np.array([1/np.prod([tval(t,i,j) for j in range(5) if j!=i]) for i in range(5)])
def A(n,i,j,k):
    l,m=[a for a in range(5) if a not in (i,j,k)]
    Nn=lambda a,b:n[eid(a,b)]; T=lambda a,b:t[eid(a,b)]
    return T(j,k)**2*Nn(j,k)*Nn(i,l)*Nn(i,m)+T(i,k)**2*Nn(i,k)*Nn(j,l)*Nn(j,m)-T(i,j)**2*Nn(i,j)*Nn(k,l)*Nn(k,m)
good=[]
for n in N:
    M=np.array([[tval(t,i,j)*n[eid(i,j)] if i!=j else 0 for j in range(5)] for i in range(5)])
    sc=np.abs(n).max()
    lin=np.abs(M@c).max()/sc
    sv=np.linalg.svd(M,compute_uv=False); rk=sv[2]/sv[0]
    D=np.prod([n[eid(0,k)]*n[eid(1,k)] for k in range(2,5)])
    sx=abs(4*t[eid(0,2)]**2*t[eid(1,2)]**2*D-A(n,0,1,2)**2-4*(t[0]*t[eid(0,2)]*t[eid(1,2)])**2)/(sc**6)
    good.append((lin,rk,sx,np.abs(n).min()/sc))
g=np.array(good); ok=(g[:,0]<1e-8)&(g[:,1]<1e-8)
print("n points",len(N),"satisfying linear+rank2:",ok.sum())
print(np.array2string(g[np.argsort(g[:,0])][:20],precision=2))
