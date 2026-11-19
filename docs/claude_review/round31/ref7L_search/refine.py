import mpmath as mp, pickle, numpy as np, itertools, sys
mp.mp.dps=60
def refine(cols,h,z0,npts=5,iters=60):
    n=len(cols)
    C=[[a[i]-a[0] for (a,e) in cols] for i in range(1,npts)]
    E=[e for a,e in cols]
    z=mp.matrix([mp.mpc(complex(x)) for x in z0])
    def F(z):
        r=[z[k] for k in range(n)]; s=[z[n+k] for k in range(n)]
        out=[]
        for m in range(1,h):
            for row in C: out.append(mp.fsum(row[k]*(r[k]**m-s[k]**m) for k in range(n)))
        out.append(mp.fsum(E[k]*(r[k]+s[k]) for k in range(n))); out.append(r[0]-s[0]-1)
        return mp.matrix(out)
    def J(z):
        r=[z[k] for k in range(n)]; s=[z[n+k] for k in range(n)]
        rows=[]
        for m in range(1,h):
            for row in C: rows.append([row[k]*m*r[k]**(m-1) for k in range(n)]+[-row[k]*m*s[k]**(m-1) for k in range(n)])
        rows.append(E+E); rows.append([1]+[0]*(n-1)+[-1]+[0]*(n-1))
        return mp.matrix(rows)
    for it in range(iters):
        f=F(z)
        if mp.norm(f)<mp.mpf(10)**(-55): break
        Jm=J(z)
        dz=mp.lu_solve(Jm.T*Jm, -(Jm.T*f)) if Jm.rows!=Jm.cols else mp.lu_solve(Jm,-f)
        z=z+dz
    return z, mp.norm(F(z))
def polys(cols,z,npts=5):
    n=len(cols); t=None
    r=[z[k] for k in range(n)]; s=[z[n+k] for k in range(n)]
    P=[]
    for i in range(npts):
        roots=[]
        for k,(a,e) in enumerate(cols):
            roots+= [r[k]]*a[i]+[s[k]]*(e-a[i])
        c=[mp.mpc(1)]
        for x in roots:
            c=[ (c[j] if j<len(c) else 0) - x*(c[j-1] if j>=1 else 0) for j in range(len(c)+1)]
        P.append(c)  # coefficients from leading
    return P
if __name__=="__main__":
    cols,Z=pickle.load(open(sys.argv[1],'rb')); h=int(sys.argv[2])
    out=[]
    for z0 in Z:
        z,res=refine(cols,h,z0)
        n=len(cols)
        roots=[z[k] for k in range(2*n)]
        sc=max(abs(x) for x in roots)
        sep=min(abs(a-b) for a,b in itertools.combinations(roots,2))/sc
        P=polys(cols,z)
        D=len(P[0])-1
        # contact: max over pairs of degree of difference = largest index j (power t^{D-j}) with nonzero coefficient diff
        maxdeg=0
        for i,j in itertools.combinations(range(5),2):
            d=[P[i][k]-P[j][k] for k in range(D+1)]
            nz=[k for k in range(D+1) if abs(d[k])>mp.mpf(10)**(-40)]
            maxdeg=max(maxdeg,D-min(nz))
        r=[z[k] for k in range(n)]; s=[z[n+k] for k in range(n)]
        realerr=max(abs(s[k]+mp.conj(r[k])) for k in range(n))
        out.append((float(sep),float(res),D-maxdeg,float(realerr),z))
        print(f"sep {float(sep):.5f} res {float(res):.1e} contact {D-maxdeg} realerr {float(realerr):.3f}")
    pickle.dump(out,open(sys.argv[1].replace('.pkl','_ref.pkl'),'wb'))
