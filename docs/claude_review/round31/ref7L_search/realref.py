import pickle,sys,numpy as np, mpmath as mp, itertools
from refine import refine, polys
cols,found=pickle.load(open(sys.argv[1],'rb')); h=int(sys.argv[2])
n=len(cols)
for sep,v in found:
    x,Y=v[:n],v[n:]
    r=(x+1j*Y)/(2j); s=(x-1j*Y)/(2j)
    z,res=refine(cols,h,np.concatenate([r,s]))
    roots=[z[k] for k in range(2*n)]
    sc=max(abs(q) for q in roots)
    sp=min(abs(a-b) for a,b in itertools.combinations(roots,2))/sc
    realerr=max(abs(z[n+k]+mp.conj(z[k])) for k in range(n))
    P=polys(cols,z); D=len(P[0])-1; md=0
    for i,j in itertools.combinations(range(5),2):
        d=[P[i][k]-P[j][k] for k in range(D+1)]
        nz=[k for k in range(D+1) if abs(d[k])>mp.mpf(10)**(-40)]
        md=max(md,D-min(nz))
    print("res",mp.nstr(res,3),"sep",mp.nstr(sp,5),"realerr",mp.nstr(realerr,3),"D",D,"contact",D-md)
    # unscaled real roots rho = 2i r
    rho=[2j*complex(z[k]) for k in range(n)]
    print("rho_k =",[complex(round(q.real,8),round(q.imag,8)) for q in rho])
    print("columns",["".join(map(str,a)) for a,e in cols])
