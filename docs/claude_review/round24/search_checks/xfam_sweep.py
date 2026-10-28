# dimension of the set swept by the family: rank of d(u_1..u_4)/d(lambda, w) along the family (log-coordinates on T_5(R))
import numpy as np
import xstruct as X
X.SIGN=-1.0
P=np.load('xfam_points.npy')
def curve_point(p,w):
    rho,r4,sig,beta,q=X.unpack(p)
    n=np.poly(np.concatenate([rho,[r4]])); In=n.imag
    qt=-q/sig; ar=In/(2*sig)
    Q2,Q1,Q0=qt; r2=np.sqrt(Q2)
    x=(Q0-w**2)/(2*r2*w-Q1); z=r2*x+w
    us=[]
    for j in (1,2,3):
        P_=np.prod([x-rho[i] for i,S in enumerate(X.S3) if j in S]); Pb=np.prod([x-np.conj(rho[i]) for i,S in enumerate(X.S3) if j in S])
        us.append(P_/Pb)
    A=np.polyval(ar,x); B=np.polyval(beta,x)
    us.append((A+1j+B*z)/(A-1j+B*z))
    return np.angle(np.array(us))   # T_5(R) = (S^1)^4 : angles
ranks=[]
for k in range(5,len(P)-5,25):
    for w in (0.37,1.3,-2.1):
        a=curve_point(P[k],w); 
        dl=(curve_point(P[k+1],w)-curve_point(P[k-1],w))
        dw=(curve_point(P[k],w+1e-6)-curve_point(P[k],w-1e-6))
        M=np.column_stack([dl/np.linalg.norm(dl),dw/np.linalg.norm(dw)])
        sv=np.linalg.svd(M,compute_uv=False); ranks.append(int(sv[-1]>1e-6))
print('rank-2 fraction:',sum(ranks),'/',len(ranks))
