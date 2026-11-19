# [J] independent check: IV.2 models are exact for an interpolating rational form Q_*; det E, disc Q_*.
from sympy import *
P,q=symbols('P q')
x=[-27,-270,135,120,40]
s=[0]*6
X=symbols('X'); poly=Poly(prod([X-xi for xi in x]),X).all_coeffs()
s1,s2,s3,s4,s5=[(-1)**k*poly[k] for k in range(1,6)]
print('s',s1,s2,s3,s4,s5)
l=97200; b0=Rational(3,2)*s5
assert l**2==-s1*s5
b=lambda i,j: x[i]*x[j]*(s1-x[i]-x[j])
h=lambda i: x[i]**3-s1*x[i]**2+s2*x[i]
B=b0*q**2; L=l*q
V=[Matrix([P*(P+q*h(i))+B*x[i], L*(P+q*h(i))]) for i in range(5)]
J={(i,j):P+q*b(i,j) for i in range(5) for j in range(i+1,5)}
Jf=lambda i,j: J[(min(i,j),max(i,j))]
n=[expand(prod([Jf(i,j) for j in range(5) if j!=i])) for i in range(5)]
for i in range(5):
  for j in range(i+1,5):
    assert expand(V[i][0]*V[j][1]-V[i][1]*V[j][0]-L*B*(x[i]-x[j])*J[(i,j)])==0
d=[expand(V[i].dot(V[i])-n[i]) for i in range(5)]
print('defects:',[Poly(di,P,q).terms()[:3] for di in d])
a,bb,c=symbols('a bb c')
Qf=lambda v: a*v[0]**2+bb*v[0]*v[1]+c*v[1]**2
sol=solve([Qf(V[i])-n[i] for i in range(3)],[a,bb,c],dict=True)[0]
sol={k:factor(v) for k,v in sol.items()}
for i in range(3,5): print('exact at',i, simplify(Qf(V[i]).subs(sol)-n[i])==0)
detE=simplify((1-sol[a])*(1-sol[c])-sol[bb]**2/4)
print('det E =',detE)
disc=factor(sol[bb]**2-4*sol[a]*sol[c]); print('disc Q* =',disc)
print('leading in P:', limit(disc*q**3/P**3,P,oo))
