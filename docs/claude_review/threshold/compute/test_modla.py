import numpy as np, time
from modla import *
rng=np.random.default_rng(1)
for p in (P1,P2):
    A=rng.integers(0,p,(37,300000)); B=rng.integers(0,p,(300000,5))
    C=mm(A,B,p)
    # exact via python ints on a few entries
    for (i,j) in [(0,0),(5,3),(36,4)]:
        v=sum(int(a)*int(b) for a,b in zip(A[i],B[:,j]))%p
        assert v==C[i,j],(v,C[i,j])
print('mm ok')
# rank test: random low-rank matrix
def naive_rank(M,p):
    M=[ [int(x)%p for x in r] for r in M]; r=0; n=len(M[0])
    for c in range(n):
        piv=None
        for i in range(r,len(M)):
            if M[i][c]: piv=i;break
        if piv is None: continue
        M[r],M[piv]=M[piv],M[r]; inv=pow(M[r][c],p-2,p); M[r]=[x*inv%p for x in M[r]]
        for i in range(len(M)):
            if i!=r and M[i][c]:
                f=M[i][c]; M[i]=[(x-f*y)%p for x,y in zip(M[i],M[r])]
        r+=1
    return r
for trial in range(5):
    p=P1
    U=rng.integers(0,p,(120,40)); V=rng.integers(0,p,(40,90))
    X=mm(U,V,p)
    X=np.vstack([X, rng.integers(0,p,(7,90))])
    rng.shuffle(X)
    E=Echelon(90,p)
    ranks=[]
    for s in range(0,X.shape[0],23):
        E.add(X[s:s+23]); ranks.append(E.rank)
    nr=naive_rank(X,p)
    assert E.rank==nr,(E.rank,nr)
    # check RREF property and kernel
    EE,pp_=E.rref(); assert np.all(EE[:,pp_]==np.eye(E.rank,dtype=np.int64))
    K=E.kernel_basis()
    assert np.all(mm(X,K.T,p)==0)
print('echelon ok', ranks)
A=rng.integers(0,P1,(3000,3000)); t=time.time(); mm(A,A,P1); print('mm 3000^3',time.time()-t)
