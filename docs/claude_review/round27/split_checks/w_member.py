# Numerical membership in W (evidence): is the image of a W-point under a torus symmetry (inversion u->1/u, or a swap
# of coordinates) again on W?  Gauss-Newton on F(a,b,1,t) ∝ P (4 equations, 3 unknowns) from random complex starts.
import numpy as np, sys
rng=np.random.default_rng(int(sys.argv[1]) if len(sys.argv)>1 else 0)
NST=int(sys.argv[2]) if len(sys.argv)>2 else 2000
def Phi(a,b,t,c=1.0):
    d=a+b+c; m=a*b*c*d; S=a*a+b*b+c*c+d*d
    za=(a*a*b+a*a*c+a*b*b+a*c*c-b*b*c-b*c*c)/2
    x=(a*(t*t+m)-2*za*t)/(t*t-m); z=za+t*(x-a)
    return np.array([(x-a)*(x-b)*(x-c)*(x+d), x**4-S/2*x*x+m+2*x*z, (x*x-b*b)*(x*x-c*c),(x*x-a*a)*(x*x-c*c),(x*x-a*a)*(x*x-b*b)])
def resid(v,P):
    F=Phi(*v); F=F/F[4]; Q=P/P[4]
    return (F[:4]-Q[:4])/(1+np.abs(Q[:4]))
def jac(v,P,h=1e-7):
    r0=resid(v,P); J=np.zeros((4,3),dtype=complex)
    for i in range(3):
        e=np.zeros(3,dtype=complex); e[i]=h; J[:,i]=(resid(v+e,P)-r0)/h
    return J
def member(P,starts):
    best=1e9
    with np.errstate(all='ignore'):
        for _ in range(starts):
            v=(rng.normal(size=3)+1j*rng.normal(size=3))*np.array([2,2,3])
            for k in range(60):
                r=resid(v,P); J=jac(v,P)
                if not (np.all(np.isfinite(r)) and np.all(np.isfinite(J))): break
                dv=np.linalg.lstsq(J,-r,rcond=None)[0]
                v=v+dv
                if not np.all(np.isfinite(v)) or np.linalg.norm(v)>1e6: break
                if np.linalg.norm(dv)<1e-13*(1+np.linalg.norm(v)): break
            if np.all(np.isfinite(v)):
                rr=np.linalg.norm(resid(v,P))
                if np.isfinite(rr) and rr<best: best=rr
            if best<1e-10: break
    return best
a0,b0,t0=0.37+0.0j,1.91+0.0j,0.53+0.0j
P=Phi(a0,b0,t0)
print("control: W point itself, residual",member(P,200),flush=True)
rnd=rng.normal(size=5)+0j
print("control: random point, residual",member(rnd,NST),flush=True)
print("inversion 1/P, residual",member(1/P,NST),flush=True)
for (i,j) in [(0,1),(0,2),(1,2)]:
    Q=P.copy(); Q[[i,j]]=Q[[j,i]]
    print(f"swap {i}<->{j}, residual",member(Q,NST),flush=True)
Q=1/P; Q[[0,1]]=Q[[1,0]]
print("inversion + swap 0<->1, residual",member(Q,NST),flush=True)
