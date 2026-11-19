# Reduced system for the (G1) polynomial rank-one family (mu=m^2 square case):
# sigma_e (10 edges of K5), c<0, K>0 ; for every vertex i:
#   prod_{e at i} sigma_e = K ,  sum_{e at i} (sigma_e + c/sigma_e) = 0
import numpy as np, itertools, sys
from scipy.optimize import least_squares
E=[(0,1),(0,2),(0,3),(0,4),(1,2),(1,3),(1,4),(2,3),(2,4),(3,4)]
inc=[[k for k,e in enumerate(E) if i in e] for i in range(5)]
rng=np.random.default_rng(int(sys.argv[1]) if len(sys.argv)>1 else 0)
c=-1.0
def res(z):
    s=z[:10]; K=z[10]
    r=[]
    for i in range(5):
        r.append(np.prod(s[inc[i]])-K)
        r.append(np.sum(s[inc[i]]+c/s[inc[i]]))
    return np.array(r)
sols=[]
for t in range(30):
    z0=np.concatenate([rng.normal(size=10)*rng.choice([0.5,1,2,4]),[abs(rng.normal())+0.1]])
    try: S=least_squares(res,z0,xtol=1e-14,ftol=1e-14,gtol=1e-14,max_nfev=300)
    except Exception: continue
    if np.linalg.norm(S.fun)<1e-10:
        s=S.x[:10]
        m_over = -S.x[10]/c**2   # (m-1)/(m+1) = c/a^2 = -K/c^2
        rho=(c-s**2)/s           # proportional to P-coordinate (times 1/2m)
        gap=np.min(np.diff(np.sort(rho)))
        sols.append((gap,S.x))
print(len(sols),'solutions; nondegenerate (min rho gap>1e-4):',sum(1 for g,_ in sols if g>1e-4))
for g,x in sorted(sols,key=lambda t:-t[0])[:10]:
    print('gap %.3g'%g, np.round(x,5))
