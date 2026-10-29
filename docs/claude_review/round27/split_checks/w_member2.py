# Numerical membership in W (evidence only), reduced to one unknown.
# W-point with c = 1: q_j = F_j/F_4. From q_2 = (x^2-1)/(x^2-a^2), q_3 = (x^2-1)/(x^2-b^2):
#   a^2 = x^2 - (x^2-1)/q_2, b^2 = x^2 - (x^2-1)/q_3.  Then q_0 = (x-1)(x+a+b+1)/((x+a)(x+b))  (one equation in x)
#   and q_1 = (x^4 - S x^2/2 + m +- 2 x z)/((x^2-a^2)(x^2-b^2)), z^2 = m x^2 + B  (check).
# A point P is (numerically) on W iff some root x of the q_0-equation (any signs of a, b) also matches q_1 (some sign of z).
import numpy as np, sys, itertools
rng=np.random.default_rng(int(sys.argv[1]) if len(sys.argv)>1 else 0)
NST=int(sys.argv[2]) if len(sys.argv)>2 else 400
def Phi(a,b,c,x,zs=+1):
    d=a+b+c; m=a*b*c*d; sq=[a*a,b*b,c*c,d*d]; S=sum(sq)
    e3=sum(sq[i]*sq[j]*sq[k] for i,j,k in itertools.combinations(range(4),3)); B=(e3-S*m)/4
    z=zs*np.sqrt(m*x*x+B+0j)
    return np.array([(x-a)*(x-b)*(x-c)*(x+d), x**4-S/2*x*x+m+2*x*z, (x*x-b*b)*(x*x-c*c),(x*x-a*a)*(x*x-c*c),(x*x-a*a)*(x*x-b*b)])
def member(P):
    q=P/P[4]; best=1e9
    with np.errstate(all='ignore'):
        for sa,sb in itertools.product([1,-1],repeat=2):
            def ab(x):
                a=sa*np.sqrt(x*x-(x*x-1)/q[2]+0j); b=sb*np.sqrt(x*x-(x*x-1)/q[3]+0j); return a,b
            def E(x):
                a,b=ab(x); return (x-1)*(x+a+b+1)-q[0]*(x+a)*(x+b)
            for _ in range(NST):
                x=(rng.normal()+1j*rng.normal())*3
                for k in range(80):
                    h=1e-7*(1+abs(x)); f=E(x); fp=(E(x+h)-f)/h
                    if not np.isfinite(f) or fp==0: break
                    dx=-f/fp; x=x+dx
                    if abs(dx)<1e-14*(1+abs(x)): break
                if not np.isfinite(x) or abs(E(x))>1e-9*(1+abs(x))**4: continue
                a,b=ab(x)
                for zs in (1,-1):
                    F=Phi(a,b,1.0,x,zs);
                    if F[4]==0: continue
                    r=np.abs(F/F[4]-q).max()/(1+np.abs(q).max())
                    best=min(best,r)
    return best
a0,b0,c0,x0=0.37,1.91,1.0,2.3
P=Phi(a0,b0,c0,x0)
print("control: W point, residual %.2e"%member(P),flush=True)
P2=Phi(1.3,-0.7,1.0,0.9+0.4j,-1)
print("control: complex W point, residual %.2e"%member(P2),flush=True)
print("control: random point, residual %.2e"%member(rng.normal(size=5)+1j*rng.normal(size=5)),flush=True)
for name,Q in [("inversion 1/P",1/P),("inversion 1/P2",1/P2)]:
    print(name,"residual %.2e"%member(Q),flush=True)
for (i,j) in [(0,1),(2,3),(0,2),(1,2)]:
    Q=P.copy(); Q[[i,j]]=Q[[j,i]]
    print(f"swap {i}<->{j}: residual %.2e"%member(Q),flush=True)
Q=1/P; Q[[0,1]]=Q[[1,0]]
print("inversion+swap 0<->1: residual %.2e"%member(Q),flush=True)
