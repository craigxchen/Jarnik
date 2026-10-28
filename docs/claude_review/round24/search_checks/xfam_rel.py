import numpy as np, itertools
import xstruct as X
P=np.load('xfam_points.npy')
def coords(p):
    rho,r4,sig,beta,q=X.unpack(p)
    return dict(q1=q[1],q0=q[2],sig=-np.exp(p[16]),b2=beta[1],b1=beta[2],b0=beta[3],
                Im1=rho[0].imag,Im2=rho[1].imag,Im3=rho[2].imag,Re1=rho[0].real,Re2=rho[1].real,Re3=rho[2].real)
C=[coords(p) for p in P]
def relation(a,b,maxdeg=8):
    A=np.array([c[a] for c in C]); B=np.array([c[b] for c in C])
    sa,sb=np.std(A),np.std(B); ma,mb=np.mean(A),np.mean(B)
    for d in range(1,maxdeg+1):
        mons=[(i,j) for i in range(d+1) for j in range(d+1-i)]
        M=np.array([[A[k]**i*B[k]**j for (i,j) in mons] for k in range(len(A))])
        M=M/np.linalg.norm(M,axis=0)
        sv=np.linalg.svd(M,compute_uv=False)
        if sv[-1]/sv[0]<1e-11: return d,sv[-1]/sv[0]
    return None,sv[-1]/sv[0]
for a,b in [('q1','sig'),('q1','q0'),('q1','b2'),('sig','q0'),('q1','Im1'),('q1','b0')]:
    print(a,b,relation(a,b))
