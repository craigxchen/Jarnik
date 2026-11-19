# [J] per-pair equivalence: deg(z_a-z_b) = D + deg D_ab - (alpha_a+alpha_b)/2, on all pairs of all 2^D points
import sympy as sp, itertools
t=sp.symbols('t'); I=sp.I
D=5; rho=[sp.Rational(r)+sp.Rational(h)*I for r,h in ((0,1),(1,2),(3,1),(-2,3),(5,1))]
cj=lambda f:sp.expand(sp.conjugate(f).subs(sp.conjugate(t),t))
pts=list(itertools.product((0,1),repeat=D))
z=lambda a:sp.expand(sp.prod([(t-p) if ak else (t-sp.conjugate(p)) for p,ak in zip(rho,a)]))
P={a:sp.expand(sp.prod([(t-p) for p,ak in zip(rho,a) if ak])) for a in pts}
bad=0;n=0
for a,b in itertools.combinations(pts,2):
    Dm=sp.expand((cj(P[a])*P[b]-P[a]*cj(P[b]))/(2*I))
    if Dm==0: continue
    lhs=sp.degree(sp.expand(z(a)-z(b)),t); rhs=D+sp.degree(Dm,t)-(2*sum(a)+2*sum(b))//2
    g=sp.gcd(sp.expand(P[a]*cj(P[a])),sp.expand(P[b]*cj(P[b])))
    n+=1; bad+= (lhs!=rhs) or sp.rem(Dm,g,t)!=0
print('pairs',n,'identity/divisibility failures',bad)
