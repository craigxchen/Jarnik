import sys, numpy as np
from scipy.optimize import least_squares
from common import *
rng = np.random.default_rng(int(sys.argv[1]) if len(sys.argv)>1 else 0)
N = int(sys.argv[2]) if len(sys.argv)>2 else 200
out = open(sys.argv[3] if len(sys.argv)>3 else 'num_lm.out','w')
stats = {'conv':0,'nondeg':0}
for it in range(N):
    L = rng.normal(size=(2,34)); c0 = rng.normal(size=2)
    def F(x):
        return np.concatenate([residual(x), L@x - c0])
    x0 = np.concatenate([rng.normal(size=15)*1.5, rng.normal(size=15)*1.5, rng.uniform(0,2*np.pi,4)])
    sol = least_squares(F, x0, method='lm', xtol=1e-15, ftol=1e-15, gtol=1e-15, max_nfev=20000)
    x = sol.x; rn = np.linalg.norm(F(x))
    if rn > 1e-9: continue
    stats['conv'] += 1
    b = x[15:30]; r = x[:15]+1j*b
    mind = min(abs(r[i]-r[k]) for i in range(15) for k in range(i+1,15))
    minc = min(abs(r[i]-np.conj(r[k])) for i in range(15) for k in range(15))
    kap = np.sin(x[30:34])
    scale = np.max(np.abs(r))
    tt = 0.37
    us = [np.polyval(Gpoly(x,j),tt)/np.conj(np.polyval(Gpoly(x,j),tt)) for j in range(4)]
    mindu = min([abs(us[i]-us[k]) for i in range(4) for k in range(i+1,4)]+[abs(u-1) for u in us])
    line = f"{it} res={rn:.1e} min|Im r|={np.min(abs(b)):.3e} minsep={mind:.3e} min|r-conj r'|={minc:.3e} min|kappa|={np.min(abs(kap)):.3e} max|r|={scale:.2e} min|u_j-u_k|={mindu:.3e}"
    nd = np.min(abs(b))>1e-4 and mind>1e-4 and np.min(abs(kap))>1e-4 and mindu>1e-6 and scale<1e6
    if nd:
        stats['nondeg'] += 1
        line += " NONDEG " + repr(list(x))
    print(line, file=out, flush=True)
print(stats, file=out)
