# (a) forced characters integral: Smith form of zero-sum row-difference lattice of K_S restricted to I.
# (b) min ||v||_1 over characters with s=||lambda||_1=4 (exhaustive, S=A or random S) vs 1.2q, and the
#     inequality (s+1)h >= q||lambda||^2 + s q/2 ; supp(v) >= (q+1)/2.
# (c) Vandermonde: K^T mu mod q has >= (q+1)/2 nonzeros for random nonzero mu supported on < q/2 indices.
import numpy as np, random, itertools, math
from sympy import Matrix, ZZ
from sympy.matrices.normalforms import smith_normal_form
def K_of(q):
    sq=set((i*i)%q for i in range(1,q)); c=np.array([0]+[1 if i in sq else -1 for i in range(1,q)])
    idx=np.arange(q); return c[(idx[None,:]-idx[:,None])%q]
rng=random.Random(5)
for q in [43,67,103]:
    K=K_of(q); A=[i for i in range(q) if K[i,0]==1]
    for trial in range(2):
        S=A if trial==0 else [i for i in A if rng.random()<0.9]
        D=[K[i,1:]-K[S[0],1:] for i in S[1:]]
        snf=smith_normal_form(Matrix(D),domain=ZZ)
        inv=[abs(snf[j,j]) for j in range(min(snf.shape)) if snf[j,j]!=0]
        print("q=%d |S|=%d rank=%d invariant factors all 1: %s"%(q,len(S),len(inv),all(x==1 for x in inv)))
    # (b) exhaustive s=4 with S=A
    KS=K[A][:,1:]; m=len(A); best=(10**9,None); minsupp=10**9; ineq_ok=True
    for a,b in itertools.combinations(range(m),2):
        for c,d in itertools.combinations([x for x in range(m) if x not in (a,b)],2):
            v=KS[a]+KS[b]-KS[c]-KS[d]; h=int(np.abs(v).sum()); minsupp=min(minsupp,int((v!=0).sum()))
            if h<best[0]: best=(h,(a,b,c,d))
            ineq_ok &= (5*h>=q*4+4*q/2)
    for a in range(m):
        for b,c in itertools.combinations([x for x in range(m) if x!=a],2):
            v=2*KS[a]-KS[b]-KS[c]; h=int(np.abs(v).sum()); minsupp=min(minsupp,int((v!=0).sum()))
            if h<best[0]: best=(h,(a,b,c))
            ineq_ok &= (5*h>=q*6+4*q/2)
    print("q=%d s=4: min ||v||_1 = %d = %.3f q (claim >=1.2q); min supp %d (claim >= %d); ineq ok %s"%(q,best[0],best[0]/q,minsupp,(q+1)//2,ineq_ok))
    # (c)
    worst=10**9
    for _ in range(3000):
        k=rng.randrange(1,(q+1)//2); sup=rng.sample(range(q),k); mu=np.zeros(q,dtype=np.int64)
        for j in sup: mu[j]=rng.randrange(1,q)
        worst=min(worst,int(((K.T@mu)%q!=0).sum()))
    print("q=%d Vandermonde: min nonzeros of K^T mu mod q = %d (claim >= %d)"%(q,worst,(q+1)//2))
