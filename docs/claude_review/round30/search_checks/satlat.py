# saturation V cap Z^D of the span of integer vectors, LLL-reduced (python-flint), and the balance test of Lemma 1.1
import flint, numpy as np, itertools
def saturation(A):
    A=[list(map(int,r)) for r in A]; D=len(A[0])
    M=flint.fmpz_mat(A); r=M.rank()
    # rational nullspace of A gives V^perp (columns), scale to integers
    X,nul=M.nullspace()
    Bp=[[int(X[i,j]) for i in range(D)] for j in range(nul)]   # rows span V^perp
    if nul==0: return [[1 if i==j else 0 for i in range(D)] for j in range(D)]
    Nbig=10**6
    rows=[[1 if i==j else 0 for i in range(D)]+[Nbig*Bp[k][j] for k in range(nul)] for j in range(D)]
    L=flint.fmpz_mat(rows).lll()
    basis=[]
    for i in range(D):
        row=[int(L[i,j]) for j in range(D+nul)]
        if all(v==0 for v in row[D:]): basis.append(row[:D])
    assert len(basis)==r, (len(basis),r)
    return basis
def balance_violation(A,L,h,cmax=2):
    B=np.array(saturation(A)); r=len(B)
    for c in itertools.product(range(-cmax,cmax+1),repeat=r):
        if not any(c): continue
        a=np.array(c)@B
        if not a.any(): continue
        if a[a>0].sum()<h or -a[a<0].sum()<h or np.abs(a).sum()<L+1: return a
    return None
