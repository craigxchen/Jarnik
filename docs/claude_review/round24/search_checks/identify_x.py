# algebraic identification of the x-structured solution (refine to 120 digits, then findpoly)
import numpy as np, ast, mpmath as mp
import xstruct as X
import verify_x as V
mp.mp.dps=130
lines=[l for l in open('xstruct_run.out') if 'NONDEG' in l]
l=lines[1]; y=np.array(ast.literal_eval(l[l.index('['):].strip().replace('np.float64','')))
y=[mp.mpf(v) for v in y]; sg=-1.0
for it in range(40):
    f=V.Fmp(y,sg); nf=mp.sqrt(sum(v*v for v in f))
    if nf<mp.mpf(10)**(-115): break
    h=mp.mpf(10)**(-60); J=mp.matrix(23,23)
    for k in range(23):
        yp=list(y); yp[k]+=h; ym=list(y); ym[k]-=h
        fp=V.Fmp(yp,sg); fm=V.Fmp(ym,sg)
        for i in range(23): J[i,k]=(fp[i]-fm[i])/(2*h)
    d=mp.lu_solve(J,mp.matrix([-v for v in f])); y=[y[k]+d[k] for k in range(23)]
print('residual',mp.nstr(nf,5))
names=['Re r1','Re r2','Re r3','Re r12','Re r13','Re r23','Re r123','Im r1','Im r2','Im r3','Im r12','Im r13','Im r23','Im r123','Re r4','Im r4','log(-sigma)','b2','b1','b0','q2','q1','q0']
vals=list(y); vals[16]=-mp.e**y[16]; names[16]='sigma'
mp.mp.dps=110
for nm,v in zip(names,vals):
    found=None
    for deg in range(1,13):
        p=mp.findpoly(v,deg,maxcoeff=10**8)
        if p: found=(deg,p); break
    print(f"{nm:12s} {mp.nstr(v,25):>32s}  minpoly:",found)
# Jacobian of the 23x23 x-structured system at the refined point
mp.mp.dps=130
h=mp.mpf(10)**(-50); J=mp.matrix(23,23)
for k in range(23):
    yp=list(y); yp[k]+=h; ym=list(y); ym[k]-=h
    fp=V.Fmp(yp,sg); fm=V.Fmp(ym,sg)
    for i in range(23): J[i,k]=(fp[i]-fm[i])/(2*h)
sv=mp.svd_r(J,compute_uv=False)
print('x-structured Jacobian singular values:',[mp.nstr(s,5) for s in sv])
