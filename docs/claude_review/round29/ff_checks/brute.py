import numpy as np, itertools, sys
from collections import Counter
def Lval(D,c): return D-1-int(np.floor(D/2+c))
def maxcluster(rho, c, tol=1e-7, ret=False):
    """rho: list of distinct upper-half-plane roots (squarefree N). returns max # transversals with equal Im p_m, m=1..L"""
    n=len(rho); D=n; L=Lval(D,c)
    if L<=0: return 2**n
    rho=np.array(rho,dtype=complex)
    A=np.array([[ (r**m).imag for r in rho] for m in range(1,L+1)])  # L x n
    # scale rows
    A=A/np.maximum(np.abs(A).max(axis=1,keepdims=True),1e-300)
    S=((np.arange(2**n)[:,None]>>np.arange(n))&1)*2-1
    V=S@A.T
    key=np.round(V/tol).astype(np.int64)
    # random projection hash
    h=key@np.random.RandomState(1).randint(1,2**20,size=L)
    cnt=Counter(h.tolist())
    best=max(cnt.values())
    if ret:
        hb=max(cnt,key=cnt.get); return best, S[h==hb]
    return best
if __name__=="__main__":
    for D in range(4,21):
        th=[np.pi*(2*k+1)/(2*D) for k in range(D)]
        rho=[np.exp(1j*t) for t in th]
        print(D, [maxcluster(rho,c) for c in (0,1)])
