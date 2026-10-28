import numpy as np, ast, mpmath as mp
import xstruct as X
from verify_x import refine, poly_from_roots, polymul
mp.mp.dps=60
lines=[l for l in open('xstruct_run.out') if 'NONDEG' in l]
l=lines[1]
y=np.array(ast.literal_eval(l[l.index('['):].strip().replace('np.float64','')))
X.SIGN=-1.0; print(np.linalg.norm(X.F(y)))
y,nf=refine(list(y),-1.0)
rho=[mp.mpc(y[k],y[7+k]) for k in range(7)]; r4=mp.mpc(y[14],y[15]); sig=-mp.e**y[16]
beta=[mp.mpf(1),y[17],y[18],y[19]]; q=[y[20],y[21],y[22]]
n=poly_from_roots(rho+[r4]); In=[mp.im(x) for x in n]
qt=[-v/sig for v in q]; ar=[v/(2*sig) for v in In]
def Hval(xv,zv,sgn):
    A=sum(ar[k]*xv**(len(ar)-1-k) for k in range(len(ar))); B=sum(beta[k]*xv**(3-k) for k in range(4))
    return A+sgn*mp.mpc(0,1)+B*zv
for i,r in enumerate(rho+[r4]):
    zz=mp.sqrt(qt[0]*r**2+qt[1]*r+qt[2])
    print(i,'H(+z)=',mp.nstr(abs(Hval(r,zz,1)),5),'H(-z)=',mp.nstr(abs(Hval(r,-zz,1)),5),' Hbar at conj:',mp.nstr(abs(Hval(mp.conj(r),mp.conj(zz),-1)),5),mp.nstr(abs(Hval(mp.conj(r),-mp.conj(zz),-1)),5))
