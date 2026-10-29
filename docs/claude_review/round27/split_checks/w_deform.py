# Deform a W-curve (two cusps at e) into curves with four smooth branches at e (numerical evidence for Prop 3.7 at x_R).
# W-curve in product form: F_j(t) ∝ lam_j prod_{places P, j in B_P}(t - t_P), passing through e at t = ±sqrt(m) with cusps.
# 4-point system: unknowns (t_P (15), lam_1..4 (lam_5 fixed), s_1..s_4); equations u_j(s_i)/u_5(s_i) = 1 (16 eqs).
# Solution set: dimension 4 + 3 (PGL_2) = 7 if the expected-dimension family exists.  Then dominance rank of
# (family, s) -> T must be 4.
import numpy as np, sympy as sp, itertools, sys
rng=np.random.default_rng(int(sys.argv[1]) if len(sys.argv)>1 else 0)
a,b,c=sp.Rational(3,2),sp.Rational(5,3),sp.Rational(7,5)
d=a+b+c; m=a*b*c*d; sq=[a*a,b*b,c*c,d*d]; S=sum(sq)
e3=sum(sq[i]*sq[j]*sq[k] for i,j,k in itertools.combinations(range(4),3)); B=(e3-S*m)/4
za=(a*a*b+a*a*c+a*b*b+a*c*c-b*b*c-b*c*c)/2
assert sp.simplify(za**2-(m*a*a+B))==0
t=sp.symbols('t')
den=t**2-m; xn=a*(t**2+m)-2*za*t; zn=za*den+t*(xn-a*den)
F=[(xn-a*den)*(xn-b*den)*(xn-c*den)*(xn+d*den), xn**4-S/2*xn**2*den**2+m*den**4+2*xn*zn*den**2,
   (xn**2-b**2*den**2)*(xn**2-c**2*den**2),(xn**2-a**2*den**2)*(xn**2-c**2*den**2),(xn**2-a**2*den**2)*(xn**2-b**2*den**2)]
polys=[sp.Poly(sp.expand(f),t) for f in F]
print("degrees",[p.degree() for p in polys])
# remove common factor (the constant place / any common root)
g=polys[0]
for p in polys[1:]: g=sp.gcd(g,p)
print("common factor degree",g.degree())
polys=[sp.Poly(sp.quo(p.as_expr(),g.as_expr()),t) for p in polys]
roots=[np.array([complex(r) for r in sp.Poly(p,t).nroots(n=30)]) for p in polys]
lead=[complex(p.LC()) for p in polys]
allr=np.concatenate(roots)
places=[]
for r in allr:
    if not any(abs(r-q)<1e-8 for q in places): places.append(r)
Bsets=[frozenset(j for j in range(5) if np.min(np.abs(roots[j]-q))<1e-8) for q in places]
print("number of places",len(places)," set sizes",sorted(len(B) for B in Bsets))
inc=np.array([[1 if j in Bp else 0 for Bp in Bsets] for j in range(5)])
lam0=np.array(lead)/lead[4]
sm=complex(sp.sqrt(m))
nP=len(places)
def unpack(w): return w[:nP], np.concatenate([w[nP:nP+4],[1.0]]), w[nP+4:nP+8]
def U(tp,lam,x): return lam*np.prod(np.where(inc==1,(x-tp)[None,:],1.0),axis=1)
def eqs(w):
    tp,lam,si=unpack(w); out=[]
    for x in si:
        u=U(tp,lam,x); out.extend(u[:4]/u[4]-1)
    return np.array(out)
def jac(w,h=1e-8):
    f0=eqs(w); J=np.zeros((len(f0),len(w)),dtype=complex)
    for k in range(len(w)):
        e=np.zeros(len(w),dtype=complex); e[k]=h; J[:,k]=(eqs(w+e)-f0)/h
    return J
# sanity: W-curve passes through e at ±sqrt(m)
w0=np.concatenate([np.array(places),lam0[:4],[sm,sm,-sm,-sm]])
print("residual of W-curve at cusp points",np.abs(eqs(w0)).max())
u=lambda x: U(np.array(places),lam0,x)
for x in (sm+1e-3,):
    print("cusp check: |u/u5-1| at distance 1e-3:",np.abs(u(x)[:4]/u(x)[4]-1).max())
# deform: split the cusps
for trial in range(6):
    eps=1e-2*(rng.normal()+1j*rng.normal())
    w=w0.copy(); w[nP+4:]+=np.array([eps,-eps,eps*1.3j,-eps*1.3j])
    w[:nP+4]+=1e-3*(rng.normal(size=nP+4)+1j*rng.normal(size=nP+4))
    with np.errstate(all='ignore'):
        for it in range(300):
            f=eqs(w)
            if np.linalg.norm(f)<1e-12: break
            dw=np.linalg.lstsq(jac(w),-f,rcond=None)[0]; w=w+dw
    f=eqs(w); tp,lam,si=unpack(w)
    sep=min(abs(x-y) for x,y in itertools.combinations(si,2))
    J=jac(w); sv=np.linalg.svd(J,compute_uv=False); rk=int(np.sum(sv>1e-7*sv[0]))
    _,_,Vh=np.linalg.svd(J); K=Vh[rk:].conj().T
    x0=0.3+0.2j
    y=lambda w,x: (lambda tp,lam,_: U(tp,lam,x)[:4]/U(tp,lam,x)[4])(*unpack(w))
    h=1e-7; cols=[(y(w+h*K[:,k],x0)-y(w,x0))/h for k in range(K.shape[1])]+[(y(w,x0+h)-y(w,x0))/h]
    dsv=np.linalg.svd(np.array(cols).T,compute_uv=False)
    print(f"trial {trial}: residual {np.linalg.norm(f):.1e}, min separation of the 4 points {sep:.2e}, solution-set dim {len(w)-rk} (expect 7), "
          f"dominance singular values {np.round(dsv/dsv[0],8)}",flush=True)
