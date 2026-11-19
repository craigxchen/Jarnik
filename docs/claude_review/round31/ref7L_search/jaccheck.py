import pickle,sys,numpy as np, mpmath as mp, itertools
from refine import refine
from newton import make
mp.mp.dps=100
cols,found=pickle.load(open(sys.argv[1],'rb')); h=int(sys.argv[2])
n=len(cols)
for sep,v in found:
    x,Y=v[:n],v[n:]
    r=(x+1j*Y)/(2j); s=(x-1j*Y)/(2j)
    z,res=refine(cols,h,np.concatenate([r,s]),iters=100)
    zc=np.array([complex(z[k]) for k in range(2*n)])
    F,J,_=make(cols,h)
    sv=np.linalg.svd(J(zc),compute_uv=False)
    print("res",mp.nstr(res,3),"J shape",J(zc).shape,"sing vals min/max",sv.min(),sv.max())
    rho=[mp.mpc(2j)*z[k] for k in range(n)]
    for k in range(n): print(k,"".join(map(str,cols[k][0])),mp.nstr(rho[k],40))
