# Independent recomputation of the four-star and five-block polynomial families of Part V:
# minors, isotropy pattern, and the DEGREE of the norm defects delta_i = |F_i|^2 - n_i.
import sympy as sp
from itertools import combinations
T=sp.symbols('T',real=True)
def det(a,b): return sp.expand(sp.im(sp.expand(sp.conjugate(a)*b)))
def nrm(a): return sp.expand(a*sp.conjugate(a))
def analyse(F,J,name):
    print('==',name)
    E=list(combinations(range(5),2))
    for (i,j) in E:
        q,r=sp.div(det(F[i],F[j]),J[(i,j)],T); assert r==0 and q.is_number and q!=0,(i,j,q,r)
    iso=[e for e in E if all(sp.rem(nrm(F[m]),J[e],T)==0 for m in e)]
    half=[(m,e) for e in E for m in e if sp.rem(nrm(F[m]),J[e],T)==0]
    print('isotropic blocks',[(a+1,b+1) for a,b in iso],' incidences',len(half))
    for i in range(5):
        n=sp.expand(sp.prod([J[e] for e in E if i in e]))
        d=sp.Poly(sp.expand(nrm(F[i])-n),T)
        print(' mark',i+1,'deg|F|^2',sp.Poly(nrm(F[i]),T).degree(),'deg n',sp.Poly(n,T).degree(),
              'deg delta', (d.degree() if not d.is_zero else '-inf (exact)'), 'LC n',sp.Poly(n,T).LC(), 'LC|F|^2',sp.Poly(nrm(F[i]),T).LC())
    return iso
# ---- four-star, generic-ish rational shape from the certificate point (v,t,k)=(3/2,2,2)
v,t,k=sp.Rational(3,2),2,2
r=t**2; f=(1-r)*(r**2-3*r+1); x=(k+sp.Rational(f)/k)/2; y=(sp.Rational(f)/k-k)/(2*t)
w=[v*(1+x)/2, v*(1-x)/2, v*(-r+y)/2, v*(-r-y)/2]
s1=sum(w); s2=sum(a*b for a,b in combinations(w,2)); s3=sum(a*b*c for a,b,c in combinations(w,3)); s4=sp.prod(w)
assert s3==s1*s2 and s1*s4!=0
G=sp.expand(sp.prod([T-sp.I*a for a in w])); A=-T/s1+sp.I
ws=[1]+w; u=[a**2 for a in ws]
F=[sp.expand(G+sp.Rational(u[i]-1)/u[i]*(u[i]+T**2)*A) for i in range(5)]
J={(i,j):T**2+u[i]*u[j] for i,j in combinations(range(5),2)}
analyse(F,J,'four-star (Part V.3), w='+str(w))
# ---- five-block family in the T chart
z=list(map(sp.Rational,("3/2","3","1/3","-32/15","-5/54")))
X=sp.symbols('x',real=True)
h={(j,kk):z[j]*z[kk] for j,kk in combinations(range(5),2)}
Gc=sp.expand(sp.prod([X+sp.I*h[(0,j)] for j in range(1,5)]))
B=(sp.Rational(4500480,1226641)*X**2+sp.Rational(1413120,1226641))+sp.I*(X**2+sp.Rational(26758,283071)*X+1)
assert det(Gc,B)==sp.expand((X-1)**6)
d=-sp.Rational(418284581,1009043070); uu=[a**2 for a in z]
L=[sp.expand(Gc+d*B*(uu[j]-uu[0])*(uu[0]+X**2/uu[j])) for j in range(5)]
F5=[sp.expand(sp.cancel(T**4*Lj.subs(X,1+1/T)/Lj.subs(X,1))) for Lj in L]
J5={e:sp.expand(((T+1)**2+h[e]**2*T**2)/(1+h[e]**2)) for e in h}
analyse(F5,J5,'five-block curve (Part V.3)')
