# continuation of the 1-parameter x-structured family (Gauss-Newton corrector, null-vector tangent), mpmath-free (float64)
import numpy as np, ast, sys
import xstruct as X
X.SIGN=-1.0
lines=[l for l in open('xstruct_run.out') if 'NONDEG' in l]
l=lines[1]; y0=np.array(ast.literal_eval(l[l.index('['):].strip().replace('np.float64','')))
def J(y,h=1e-7):
    return np.array([(X.F(y+h*e)-X.F(y-h*e))/(2*h) for e in np.eye(len(y))]).T
def correct(y):
    for it in range(30):
        f=X.F(y); Jm=J(y); d=np.linalg.lstsq(Jm,-f,rcond=None)[0]; y=y+d
        if np.linalg.norm(d)<1e-13: break
    return y,np.linalg.norm(X.F(y))
y0,r=correct(y0)
pts=[y0.copy()]
for direction in (1,-1):
    y=y0.copy(); t_prev=None; h=0.02
    for step in range(int(sys.argv[1]) if len(sys.argv)>1 else 300):
        _,_,Vt=np.linalg.svd(J(y)); t=Vt[-1]
        if t_prev is None: t=t*direction
        elif t@t_prev<0: t=-t
        yn,r=correct(y+h*t)
        if r>1e-10 or not np.all(np.isfinite(yn)):
            h/=2
            if h<1e-6: break
            continue
        rho,r4,sig,beta,q=X.unpack(yn)
        pts_all=np.concatenate([rho,[r4]]); allv=np.concatenate([pts_all,np.conj(pts_all)])
        iu=np.triu_indices(16,1); sep=np.min(np.abs(allv[:,None]-allv[None,:])[iu])
        if sep<1e-4 or np.max(np.abs(yn))>1e4: print('stop: degenerate/large at step',step,'sep',sep); break
        t_prev=t; y=yn; pts.append(y.copy()); h=min(h*1.3,0.05)
np.save('xfam_points.npy',np.array(pts))
P=np.array(pts); rho,r4,sig,beta,q=X.unpack(P[0])
print('points',len(P))
print('range sigma:',(-np.exp(P[:,16])).min(),(-np.exp(P[:,16])).max(),' range q1:',P[:,21].min(),P[:,21].max(),' q0:',P[:,22].min(),P[:,22].max())
print('Im n4 along family (max abs):',max(abs(np.poly(np.concatenate([X.unpack(p)[0],[X.unpack(p)[1]]]))[4].imag) for p in P))
print('q2 along family:',P[:,20].min(),P[:,20].max())
