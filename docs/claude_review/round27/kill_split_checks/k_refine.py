# Find a 4-branch solution at a given class (numpy), refine to 60 digits with mpmath on an augmented square system,
# then compute the dominance matrix (rank 4 test) at high precision. Numerical (not interval-certified).
import numpy as np, itertools, sys, time, mpmath as mp
_argv=sys.argv; sys.argv=['x','0']
exec(open('k_cover.py').read().split("classes=onesided_classes()")[0])
sys.argv=_argv
sig=eval(sys.argv[1]); seed=int(sys.argv[2]); rng=np.random.default_rng(seed)
cls=[F for F in onesided_classes() if sorted(len(B) for B in F)==sig]
F,A,D=setup(cls[int(sys.argv[3]) if len(sys.argv)>3 else 0])
print("class F:",[sorted(B) for B in F])
while True:
    s4=complex(*rng.normal(size=2))*1.5
    z=np.concatenate([rng.normal(size=15)+1j*rng.normal(size=15), rng.normal(size=4)+1j*rng.normal(size=4)])*1.5
    z,nf,it=newton(z,s4,A,D)
    if nf<1e-11:
        r=z[:15]; S=np.array([0,1,-1,s4])
        sep=min(min(abs(r[i]-r[j]) for i in range(15) for j in range(i)), min(abs(si-rj) for si in S for rj in r))
        if sep>1e-2 and min(abs(z[15:19]))>1e-2 and min(abs(z[15:19]-1))>1e-2: break
mp.mp.dps=60
X0=[mp.mpc(complex(v)) for v in list(z)+[s4]]
L=[[mp.mpc(complex(*rng.normal(size=2))) for _ in range(20)] for _ in range(4)]
Lc=[mp.fsum(L[k][j]*X0[j] for j in range(20)) for k in range(4)]
def G(X):
    r=X[:15]; c=X[15:19]; Sv=[mp.mpc(0),mp.mpc(1),mp.mpc(-1),X[19]]
    out=[]
    for a in range(4):
        for si in Sv:
            out.append(c[a]*mp.fprod([si-r[j] for j in A[a]])-mp.fprod([si-r[j] for j in D[a]]))
    for k in range(4): out.append(mp.fsum(L[k][j]*X[j] for j in range(20))-Lc[k])
    return mp.matrix(out)
def JG(X):
    h=mp.mpf('1e-30'); J=mp.matrix(20,20)
    for j in range(20):
        Xp=list(X); Xm=list(X); Xp[j]+=h; Xm[j]-=h
        col=(G(Xp)-G(Xm))/(2*h)
        for i in range(20): J[i,j]=col[i]
    return J
X=X0
for itn in range(8):
    g=G(X); ng=mp.norm(g); print("  newton",itn,"|G| =",mp.nstr(ng,5))
    if ng<mp.mpf('1e-55'): break
    dX=mp.lu_solve(JG(X),-g); X=[X[j]+dX[j] for j in range(20)]
J=JG(X); sv=mp.svd_c(J,compute_uv=False); print("  smallest singular values of augmented Jacobian:",[mp.nstr(v,4) for v in sorted(sv)[:3]])
# tangent space of the 4-dim solution variety: kernel of the 16x20 equation Jacobian
Jinv=mp.inverse(J)
K=[[Jinv[j,16+k] for j in range(20)] for k in range(4)]
J16=mp.matrix(16,20)
for i in range(16):
    for j in range(20): J16[i,j]=J[i,j]
print("  check J16*K ~ 0:",mp.nstr(max(abs(mp.fsum(J16[i,j]*K[k][j] for j in range(20))) for i in range(16) for k in range(4)),3))
print("  16x20 singular values (smallest 2):",[mp.nstr(v,4) for v in sorted(mp.svd_c(J16,compute_uv=False))[:2]])
def Phi(X,s):
    r=X[:15]; c=list(X[15:19])+[mp.mpc(1)]
    u=[c[a]*mp.fprod([s-r[j] for j,B in enumerate(F) if a in B]) for a in range(5)]
    return [mp.log(u[a]/u[4]) for a in range(4)]
s=mp.mpc('0.37','0.21'); h=mp.mpf('1e-28'); cols=[]
for k in range(4):
    Xp=[X[j]+h*K[k][j] for j in range(20)]; Xm=[X[j]-h*K[k][j] for j in range(20)]
    cols.append([(p-m)/(2*h) for p,m in zip(Phi(Xp,s),Phi(Xm,s))])
cols.append([(p-m)/(2*h) for p,m in zip(Phi(X,s+h),Phi(X,s-h))])
M=mp.matrix(4,5)
for i in range(4):
    for j in range(5): M[i,j]=cols[j][i]
print("  dominance singular values (60 digits):",[mp.nstr(v,6) for v in mp.svd_c(M,compute_uv=False)])
r=X[:15]; print("  min root separation:",mp.nstr(min(abs(r[i]-r[j]) for i in range(15) for j in range(i)),4),"s4=",mp.nstr(X[19],8))
