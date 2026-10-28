# Independent check of Prop 1.5 (lower bound) and Lemma 1.4 (upper bound) for the all-ideal Hecke large sieve.
# Lambda*(X,K) = largest eigenvalue of the (2K+1)x(2K+1) Toeplitz matrix T(k-l)=sum_{N alpha<=X} e^{4i(k-l) vartheta(alpha)}.
import numpy as np, math, sys
def ideal_angles(X):
    r=math.isqrt(X); out=[]
    for x in range(1,r+1):
        ymax=math.isqrt(max(X-x*x,0))
        y=np.arange(0,ymax+1)
        out.append(np.arctan2(y,x))
    return np.concatenate(out)
for X in [100,1000,10000,100000]:
    th=ideal_angles(X); nI=len(th)
    for K in [0,1,5,20,100]:
        m=np.arange(-2*K,2*K+1)
        Tm=np.array([np.exp(4j*k*th).sum() for k in m])
        T=np.array([[Tm[(k-l)+2*K] for l in range(2*K+1)] for k in range(2*K+1)])
        lam=np.linalg.eigvalsh(T).max()
        lb=max(nI,(2*K+1)*math.isqrt(X))
        ub=56*(X+K*X**0.5) if K>=1 else float('nan')
        print(f"X={X:7d} K={K:4d} #I={nI:7d} Lambda*={lam:12.1f} lb={lb:9d} lb_ok={lam>=lb-1e-6} ub_ok={K==0 or lam<=ub} ratio_to_(X+K sqrtX)={lam/(X+K*X**0.5):.3f}")
