import pickle,numpy as np, mpmath as mp
from refine import refine
mp.mp.dps=400
cols,found=pickle.load(open('real_11_6_1.pkl','rb'))
n=len(cols); sep,v=found[0]
x,Y=v[:n],v[n:]
z,res=refine(cols,6,np.concatenate([(x+1j*Y)/(2j),(x-1j*Y)/(2j)]),iters=200)
print("res",mp.nstr(res,3))
rho0=mp.mpc(2j)*z[0]; rho1=mp.mpc(2j)*z[1]
# gauge-invariant real quantity: ratio (x1-x0)/Y0 with Y0=1 normalisation
q=mp.re((rho1-rho0)/mp.im(rho0))
for d in [2,4,6,8,10,12,16,20,24,30]:
    p=mp.findpoly(q,d,maxcoeff=10**12)
    print(d,p)
    if p: break
