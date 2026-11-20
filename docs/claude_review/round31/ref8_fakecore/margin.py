# Pair margin w(K_i-K_i') - W/2 - 2 log P over all pairs of A, for q-1 consecutive split primes >= P=q^d,
# randomly permuted.  Reports sigma (log-width of the prime interval) and min/mean of the margin excess.
import numpy as np, math, sys, random
from sympy import nextprime
def K_of(q):
    sq=set((i*i)%q for i in range(1,q)); c=np.array([0]+[1 if i in sq else -1 for i in range(1,q)])
    idx=np.arange(q); return c[(idx[None,:]-idx[:,None])%q]
for q,d in [(383,1.99),(1019,1.99),(2003,1.99),(1019,1.7)]:
    P=int(q**d); ps=[]; p=P
    while len(ps)<q-1:
        p=nextprime(p)
        if p%4==1: ps.append(p)
    logs=np.log(np.array(ps,dtype=float)); sigma=logs[-1]-logs[0]
    K=K_of(q); A=[i for i in range(q) if K[i,0]==1]; KA=K[A][:,1:]
    rng=np.random.default_rng(7); worst=[]
    for t in range(3):
        lp=rng.permutation(logs)
        Wh=lp.sum()
        # w(K_i-K_i') = sum |K_ij-K_i'j| log p_j
        Dm=lp.sum()-(KA*lp)@KA.T.astype(float)  # |x-y|=1-xy off the diagonal
        iu=np.triu_indices(len(A),1); ex=Dm[iu]-Wh-2*math.log(P)
        worst.append((ex.min(),ex.mean()))
    print("q=%d d=%.2f sigma=%.4f sigma*sqrt(q log q)=%.3f  (min,mean) of w-W/2-2logP over pairs:"%(q,d,sigma,sigma*math.sqrt(q*math.log(q))),[(round(a,3),round(b,3)) for a,b in worst])
