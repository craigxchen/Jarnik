# (G1) polynomial ansatz: V_i = (x_i P + y_i, L x_i), x_i = u_i P + v_i, y_i const.
# Want |V_i|^2 - mu*y_i^2 = n_i = prod_{j!=i} J_ij, with L*(y_i x_j - y_j x_i) = t_ij J_ij
# Equivalent: (x_iP+y_i)^2 + L2 x_i^2 - mu y_i^2 = c_i * prod_{j!=i}(y_i x_j - y_j x_i)  (c_i free const)
# 25 polynomial-coefficient equations. Numerical least squares from random starts.
import numpy as np, sys
from numpy.polynomial import polynomial as Pn
from scipy.optimize import least_squares
rng=np.random.default_rng(int(sys.argv[1]) if len(sys.argv)>1 else 0)
def unpack(z):
    u=np.concatenate([[1.0],z[0:4]]); v=z[4:9]; y=np.concatenate([[1.0],z[9:13]])
    L2=z[13]; mu=z[14]; c=z[15:20]; return u,v,y,L2,mu,c
def res(z):
    u,v,y,L2,mu,c=unpack(z); out=[]
    for i in range(5):
        xi=np.array([v[i],u[i]])
        X=Pn.polyadd(Pn.polymul(xi,[0,1]),[y[i]])
        lhs=Pn.polyadd(Pn.polymul(X,X),Pn.polyadd(L2*Pn.polymul(xi,xi),[-mu*y[i]**2]))
        r=np.array([1.0])
        for j in range(5):
            if j==i: continue
            d=np.array([y[i]*v[j]-y[j]*v[i], y[i]*u[j]-y[j]*u[i]])
            r=Pn.polymul(r,d)
        r=c[i]*r
        lhs=np.concatenate([lhs,np.zeros(5-len(lhs))]); r=np.concatenate([r,np.zeros(5-len(r))])
        out.extend(lhs-r)
    return np.array(out)
def degen(z):
    u,v,y,L2,mu,c=unpack(z)
    lead=[];roots=[]
    for i in range(5):
        for j in range(i+1,5):
            a=y[i]*u[j]-y[j]*u[i]; b=y[i]*v[j]-y[j]*v[i]
            lead.append(abs(a)); roots.append(b/a if abs(a)>1e-12 else 1e12)
    roots=np.sort(roots)
    return min(lead), np.min(np.diff(roots)), min(abs(y)), abs(mu)
best=[]
for trial in range(400):
    z0=rng.normal(size=20)*2
    try:
        s=least_squares(res,z0,xtol=1e-15,ftol=1e-15,gtol=1e-15,max_nfev=4000)
    except Exception as e: continue
    r=np.linalg.norm(s.fun)
    d=degen(s.x)
    if r<1e-8: best.append((r,d,s.x))
print('converged to zero:',len(best))
for r,d,x in best[:40]:
    print('%.1e'%r,' minlead %.2e  minrootgap %.2e  min|y| %.2e |mu| %.2e'%d)
