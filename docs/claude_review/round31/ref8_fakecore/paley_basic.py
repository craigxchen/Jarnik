# Exact checks of Paley matrix facts and the skew conference matrix C.
import numpy as np, sys
def chi_table(q):
    sq=set((i*i)%q for i in range(1,q))
    return np.array([0]+[1 if i in sq else -1 for i in range(1,q)],dtype=np.int64)
def K_of(q):
    c=chi_table(q); idx=np.arange(q)
    return c[(idx[None,:]-idx[:,None])%q]
for q in [43,67,127,131,383]:
    K=K_of(q); I=np.eye(q,dtype=np.int64); J=np.ones((q,q),dtype=np.int64)
    assert (K.T==-K).all() and (K@K.T==q*I-J).all()
    assert (K.sum(0)==0).all() and (K.sum(1)==0).all()
    A=[i for i in range(q) if K[i,0]==1]; assert len(A)==(q-1)//2
    # row distances full rows
    D=np.abs(K[:,None,:]-K[None,:,:]).sum(2); off=~np.eye(q,dtype=bool)
    assert (D[off]==q+1).all()
    # coefficients |K_ij-K_i'j|-1 on I sum to 2 for i,i' in A
    KA=K[A][:,1:]; s=(np.abs(KA[:,None,:]-KA[None,:,:])-1).sum(2)
    offA=~np.eye(len(A),dtype=bool); assert (s[offA]==2).all()
    # C
    n=q+1; C=np.zeros((n,n),dtype=np.int64); C[:q,:q]=K; C[:q,q]=1; C[q,:q]=-1
    assert (C@C==-q*np.eye(n,dtype=np.int64)).all()
    # frame identity on I restricted to A: sum_{j in I} K_kj K_ij = q delta - 2
    G=KA@KA.T; assert (G==q*np.eye(len(A),dtype=np.int64)-2).all()
    print("q",q,"ok: K^T=-K, KK^T=qI-J, row dist q+1, |A|=(q-1)/2, margin coeff sum 2, C^2=-qI, frame identity on I")
