# Krawczyk certification for base-form zero-sum systems.
# Unknowns z = (x_1..x_D, y_1..y_D), q_k = x_k + i y_k.  Equations: for each set F in Fs (independent),
# sum_{k in F} Im q_k^m = 0, m=1..L; gauge sum x = 0; sum(x^2+y^2) = D; plus fixed coords z_j = value (rational).
# Krawczyk: K(X) = z0 - Y F(z0) + (I - Y J(X))(X - z0) inside interior X  => unique zero in X.
import mpmath as mp
from mpmath import iv
import numpy as np
mp.mp.dps=60; iv.dps=60
def powparts(x,y,L,ctx):
    # returns lists re[m], im[m] of (x+iy)^m and derivative helpers
    re=[ctx.mpf(1)]; im=[ctx.mpf(0)]
    for m in range(1,L+1):
        re.append(re[-1]*x-im[-1]*y); im.append(re[-1-0]*0+im[-1]*x+re[-2+1]*y if False else im[-1]*x+re[-2]*y if False else None)
    return re,im
def pw(x,y,L):
    re=[1]; im=[0]
    for m in range(1,L+1):
        a=re[-1]*x-im[-1]*y; b=re[-1]*y+im[-1]*x
        re.append(a); im.append(b)
    return re,im
def Ffun(z,Fs,L,D,fix,one):
    x=z[:D]; y=z[D:]
    P=[pw(x[k],y[k],L) for k in range(D)]
    out=[]
    for F in Fs:
        for m in range(1,L+1):
            s=0*one
            for k in range(D):
                if F[k]: s=s+P[k][1][m]
            out.append(s)
    out.append(sum(x[k] for k in range(D)))
    out.append(sum(x[k]**2+y[k]**2 for k in range(D))-D)
    for (j,v) in fix: out.append(z[j]-v)
    return out
def Jfun(z,Fs,L,D,fix,one):
    # d Im q^m / dx = m Im q^{m-1}, / dy = m Re q^{m-1}
    x=z[:D]; y=z[D:]; n=2*D
    P=[pw(x[k],y[k],L) for k in range(D)]
    rows=[]
    for F in Fs:
        for m in range(1,L+1):
            row=[0*one]*n
            for k in range(D):
                if F[k]:
                    row[k]=m*P[k][1][m-1]; row[D+k]=m*P[k][0][m-1]
            rows.append(row)
    rows.append([one]*D+[0*one]*D)
    rows.append([2*x[k] for k in range(D)]+[2*y[k] for k in range(D)])
    for (j,v) in fix:
        row=[0*one]*n; row[j]=one; rows.append(row)
    return rows
def refine(z0,Fs,L,fix,iters=60):
    D=len(z0)//2; z=mp.matrix([mp.mpf(float(t)) for t in z0])
    for it in range(iters):
        F=mp.matrix(Ffun(list(z),Fs,L,D,fix,mp.mpf(1)))
        J=mp.matrix(Jfun(list(z),Fs,L,D,fix,mp.mpf(1)))
        dz=mp.lu_solve(J,F); z=z-dz
        if mp.norm(dz)<mp.mpf(10)**(-50): break
    return z, mp.norm(F)
def krawczyk(z,Fs,L,fix,rad=mp.mpf(10)**-30):
    D=len(z)//2; n=2*D
    J0=mp.matrix(Jfun(list(z),Fs,L,D,fix,mp.mpf(1)))
    Y=mp.inverse(J0)
    X=[iv.mpf([z[i]-rad,z[i]+rad]) for i in range(n)]
    Fz=Ffun([iv.mpf(z[i]) for i in range(n)],Fs,L,D,[(j,iv.mpf(v)) for j,v in fix],iv.mpf(1))
    JX=Jfun(X,Fs,L,D,[(j,iv.mpf(v)) for j,v in fix],iv.mpf(1))
    ok=True
    for i in range(n):
        s=iv.mpf(z[i])
        for j in range(n): s=s-iv.mpf(Y[i,j])*Fz[j]
        for j in range(n):
            a=iv.mpf(1 if i==j else 0)
            for l in range(n): a=a-iv.mpf(Y[i,l])*JX[l][j]
            s=s+a*(X[j]-iv.mpf(z[j]))
        if not (s.a>X[i].a and s.b<X[i].b): ok=False; break
    return ok
def certify(q,Fs,L,nfix=None):
    """Fs: independent zero-sum sets (list of 0/1 lists). Chooses fixed coordinates to make the system square."""
    D=len(q); n=2*D; z0=np.concatenate([q.real,q.imag])
    neq=len(Fs)*L+2; k=n-neq
    if k<0: raise ValueError("overdetermined: use independent sets only / not square")
    # choose fixed coordinates greedily to make the Jacobian well conditioned
    import itertools
    Jn=np.array([[float(t) for t in row] for row in Jfun(list(z0),Fs,L,D,[],1.0)])
    fixidx=[]
    for _ in range(k):
        best=None
        for j in range(n):
            if j in fixidx: continue
            e=np.zeros(n); e[j]=1
            M=np.vstack([Jn]+[np.eye(n)[i] for i in fixidx]+[e])
            sv=np.linalg.svd(M,compute_uv=False)[-1] if M.shape[0]>=n else np.linalg.matrix_rank(M)
            val=np.linalg.matrix_rank(M,tol=1e-8)
            if best is None or val>best[0]: best=(val,j)
        fixidx.append(best[1])
    fix=[(j,mp.mpf(round(z0[j]*1000))/1000) for j in fixidx]   # rational fixed values
    z,res=refine(z0,Fs,L,fix)
    ok=krawczyk(z,Fs,L,fix)
    qq=[complex(z[i],z[D+i]) for i in range(D)]
    return ok,res,qq,fix,z
