# Sanity check of (II.2) and the Gram bound on fake binary profiles: rows of a Sylvester Hadamard
# matrix (+-1 layers of equal weight; pair distance exactly W/2) and of random sign rows with d_ij>W/2.
import numpy as np, itertools, random
rng=np.random.default_rng(1)
def syl(k):
    H=np.array([[1]])
    for _ in range(k): H=np.block([[H,H],[H,-H]])
    return H
for k in (5,6,7,8):
    S=syl(k)[:, :]                  # m rows, layers = columns, equal weights
    S=S[:,1:]                       # drop constant column: pair distances = 2^{k-1} > W/2=(2^k-1)/2
    m,Lw=S.shape
    V=S/np.sqrt(Lw)
    G=V@V.T
    off=G-np.eye(m)
    assert off.max()<=1e-12
    nrm=np.linalg.norm(V,2)**2
    # sample ordered 5-tuples
    T=rng.integers(0,m,size=(200000,5))
    sub=S[T]                        # (n,5,Lw)
    const=(np.abs(sub.sum(axis=1))==5).mean(axis=1)
    frac=1-const
    print(f"m={m}: ||V||^2={nrm:.3f} (<=2), mean logH_S/W={frac.mean():.5f}, P(|.-15/16|>0.02)={np.mean(np.abs(frac-15/16)>0.02):.4f}")
