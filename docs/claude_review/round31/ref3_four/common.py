import numpy as np, itertools
E=[(i,j) for i in range(5) for j in range(i+1,5)]
EI={e:k for k,e in enumerate(E)}
def eid(i,j): return EI[(min(i,j),max(i,j))]
def tval(t,i,j):  # antisymmetric extension
    return t[eid(i,j)] if i<j else -t[eid(i,j)]
# log-linear map: log D_ij = sum_{k!=i,j} (x_ik + x_jk)
Madj=np.zeros((10,10))
for a,(i,j) in enumerate(E):
    for k in range(5):
        if k in (i,j): continue
        Madj[a,eid(i,k)]+=1; Madj[a,eid(j,k)]+=1
def real_sample(rng,t=None,theta=None):
    if theta is None:
        theta=np.sort(rng.uniform(0,np.pi,5)); theta-=theta[4]  # rotation: theta5=0
    if t is None:
        t=np.array([np.sign(np.sin(theta[j]-theta[i]))*rng.integers(1,6) for (i,j) in E],float)
    D=np.array([t[a]**2/np.sin(theta[j]-theta[i])**2 for a,(i,j) in enumerate(E)])
    x=np.linalg.solve(Madj,np.log(D)); n=np.exp(x)
    S=np.array([np.prod([n[eid(i,j)] for j in range(5) if j!=i]) for i in range(5)])
    P=np.sqrt(S)*np.exp(1j*theta)
    return t,n,P,theta
def Svec(n): return np.array([np.prod([n[eid(i,j)] for j in range(5) if j!=i]) for i in range(5)])
def det2(a,b): return (np.conj(a)*b).imag
