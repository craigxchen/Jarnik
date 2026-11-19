"""Exact polynomial analysis for a balanced five-mark cubic pencil with a short cycle row."""
from itertools import combinations
from math import lcm,prod
import sympy as s
Y=(1,2,-9,13,-7)
x=[s.Rational(1,y) for y in Y]
S=sum(x); e2=sum(a*b for a,b in combinations(x,2))
z=[v**3-S*v**2 for v in x]
E=list(combinations(range(5),2)); F=[(i,(i+1)%5) for i in range(5)]
b=lambda i,j:s.cancel((x[i]*z[j]-x[j]*z[i])/(x[j]-x[i]))
ell=lambda i,j:s.cancel((z[j]-z[i])/(x[j]-x[i]))
P=s.symbols('P'); J=lambda i,j:P+b(i,j)
ds=[b(i,j) for i,j in E]; assert len(set(ds))==10
D=prod(J(i,j) for i,j in F)
A=[prod(J(i,j) for i,j in F if k in (i,j)) for k in range(5)]
N=[s.cancel(D*prod(J(i,j) for i,j in E if k in (i,j))/A[k]**2) for k in range(5)]
def rows(f,g):
 return [f.subs(P,-b(i,j))+ell(i,j)*g.subs(P,-b(i,j)) for i,j in F]
hs=s.symbols('h:3');is_=s.symbols('i:2'); f=sum(hs[i]*P**i for i in range(3));g=sum(is_[i]*P**i for i in range(2))
ss=list(s.linsolve(rows(f,g)+[hs[2]-1],hs+is_));print('short',ss)
if not ss:quit()
h=f.subs(dict(zip(hs+is_,ss[0]))); I=g.subs(dict(zip(hs+is_,ss[0])))
r=[s.cancel((h*v+I*(-P+u))/a) for v,u,a in zip(x,z,A)]; print('r',r)
# T_i numerator F x_i+G(-P+z_i) with G leading -P^3.
# leading sqrt if alpha norm=P+c: 1/2 (sum_F b+star_sum b−c) = const−z_i−e2 xi;
# hence F3=−e2. G2 and c coupled.
fp=s.symbols('f:3');gp=s.symbols('g:3')
f=-e2*P**3+sum(fp[i]*P**i for i in range(3))
g=-P**3+sum(gp[i]*P**i for i in range(3))
ss=list(s.linsolve(rows(f,g),fp+gp));print('long solution',ss)
if not ss:quit()
ss=ss[0];subs=dict(zip(fp+gp,ss));u=f.subs(subs).subs(gp[0],0).subs(gp[1],0).subs(gp[2],0)
v=g.subs(subs).subs(gp[0],0).subs(gp[1],0).subs(gp[2],0)
# Replace all remaining free variables with zero consistently.
symbols=(u.free_symbols|v.free_symbols)-{P};simple={q:0 for q in symbols};u=u.subs(simple);v=v.subs(simple)
T=[s.cancel((u*a+v*(-P+zz))/aa) for a,zz,aa in zip(x,z,A)]
kappa=sum(b(i,j) for i,j in F)+2*s.expand(v).coeff(P,2)
print('c',kappa)
Bk=[s.Poly(s.expand((P+kappa)*vv**2-n),P) for vv,n in zip(T,N)]
print('orders',[bb.degree() for bb in Bk])
la=[bb.nth(3) for bb in Bk];print('lambda/r',list(set(aa/rr for aa,rr in zip(la,r) if rr)))
print('quadratic interpolation',s.interpolate(list(zip(x,la)),P).factor())
# Fixed rational normalization produces determinant exactly D.
assert s.cancel((u*I-v*h)/D)==1
assert kappa==s.Rational(47,234)
assert all(s.cancel(rr-a-s.Rational(13,7))==0 for a,rr in zip(x,r))
tau=s.Rational(2677,2*415233)
T=[vv+tau*rr for vv,rr in zip(T,r)]
assert all(s.Poly(s.cancel((P+kappa)*vv**2-n),P).degree()==2 for vv,n in zip(T,N))
L=lcm(*(int(y.q) for y in r+x))
ri=[int(rr*L) for rr in r]
for i,j in E:
 assert s.cancel(T[i]*ri[j]-T[j]*ri[i]-L*(x[j]-x[i])*D*J(i,j)/(A[i]*A[j]))==0
print('integral r multiplier',L,'r',ri)
