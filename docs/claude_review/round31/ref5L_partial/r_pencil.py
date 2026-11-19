# Independent check of the IV.2 pencil isotropy failure (Part V.3, first claim)
import sympy as sp
from itertools import combinations
P,q=sp.symbols('P q')
x=[-27,-270,135,120,40]
X=sp.symbols('X')
poly=sp.expand(sp.prod([X-xi for xi in x]))
c=sp.Poly(poly,X).all_coeffs()  # X^5 - s1 X^4 + s2 X^3 - s3 X^2 + s4 X - s5
s1,s2,s3,s4,s5=-c[1],c[2],-c[3],c[4],-c[5]
print('s',s1,s2,s3,s4,s5)
l=sp.sqrt(-s1*s5); b0=sp.Rational(3,2)*s5
print('l',l,'b0',b0)
assert l.is_Integer and b0.is_Integer and s4==0
b={}
for i,j in combinations(range(5),2): b[i,j]=b[j,i]=x[i]*x[j]*(s1-x[i]-x[j])
h=[xi**3-s1*xi**2+s2*xi for xi in x]
assert len({b[e] for e in combinations(range(5),2)})==10
B=b0*q**2; L=l*q
V=[(P*(P+q*h[i])+B*x[i], L*(P+q*h[i])) for i in range(5)]
for i,j in combinations(range(5),2):
    det=sp.expand(V[i][0]*V[j][1]-V[i][1]*V[j][0])
    assert sp.expand(det-L*B*(x[i]-x[j])*(P+q*b[i,j]))==0
    for k in (i,j):
        nrm=sp.expand(V[k][0]**2+V[k][1]**2)
        red=sp.expand(nrm.subs(P,-q*b[i,j]))
        br=(-b[i,j]*(h[k]-b[i,j])+b0*x[k])**2+l**2*(h[k]-b[i,j])**2
        assert sp.expand(red-q**4*br)==0 and br>0
print('symbolic remainder identity OK on all 20 incidences; brackets positive')
# integer instance with P,q, blocks need not be prime for this test
import math
qq=2*3*5*7*11*13*31*61*389*1481  # include primes of l,b0,x
for pr in sp.factorint(int(l*b0*sp.prod(x))): assert qq%pr==0, pr
PP=10**60+7
iso=0
for i,j in combinations(range(5),2):
    J=PP+qq*b[i,j]
    for k in (i,j):
        a0=PP*(PP+qq*h[k])+b0*qq**2*x[k]; a1=l*qq*(PP+qq*h[k])
        r=(a0*a0+a1*a1)%J
        iso+= (r==0)
        # also prime-level: gcd
        g=math.gcd(int(a0*a0+a1*a1),int(J))
print('integer instance: isotropic whole-block incidences =',iso)
