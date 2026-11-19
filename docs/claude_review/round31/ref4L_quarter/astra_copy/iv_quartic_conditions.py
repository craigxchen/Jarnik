"""Exact shell discrepancy after a cubic star interpolation."""
from sympy import symbols, expand, rem, Poly, factor
x,y=symbols("x y")
s,t,u,v,w=symbols("s t u v w")
G=x**5-s*x**4+t*x**3-u*x*x+v*x-w
h=x**3-s*x*x+t*x
co=[1,-s,t,-u,v,-w]
power=[5]
for k in range(1,10):
    power.append(expand(-sum(co[i]*power[k-i] for i in range(1,min(k,6))) -(k*co[k] if k<6 else 0)))
b=x*y*(s-x-y)
def mod(a): return rem(expand(a),G,x)
def trace(a):
    p=Poly(a,y)
    return mod(sum(v*power[k[0]] for k,v in p.terms())-a.subs(y,x))
q=[4]+[trace(b**k) for k in range(1,5)]
es=[1]
for k in range(1,5):
    es.append(mod(sum((-1)**(r-1)*es[k-r]*q[r] for r in range(1,k+1))/k))
assert mod(es[1]-2*h)==0 and mod(h*h-es[2]-v*x*x+(s*v+3*w)*x-s*w)==0
R1=mod(3*h*w*x-2*s*w*h-es[3]).subs(v,0)
R0=mod(9*w*w*x*x/4-s*w*h*h-es[4]).subs(v,0)
print('linear-in-p defect',factor(R1))
print('constant-in-p defect',factor(R0))

from itertools import combinations
from math import prod,isqrt,gcd
xx=(-27,-270,135,120,40)
elem=[1]+[sum(prod(a) for a in combinations(xx,j)) for j in range(1,6)]
S,E2,E3,E4,E5=elem[1:]
assert E4==0 and -S*E5==isqrt(-S*E5)**2
baseL=isqrt(-S*E5)
assert E5%2==0
baseB=3*E5//2
edges=list(combinations(range(5),2))
bs={e:xx[e[0]]*xx[e[1]]*(S-xx[e[0]]-xx[e[1]]) for e in edges}
assert len(set(bs.values()))==10
h=[z**3-S*z*z+E2*z for z in xx]
print('elementary symmetrics',elem,'L/A',baseL,'B/A²',baseB)
print('pair offsets', list(bs.values()))
p,A=symbols('p A')
Y=[baseL*A*(p+A*hj) for hj in h]
X=[p*(p+A*hj)+baseB*A**2*z for hj,z in zip(h,xx)]
for (i,j),bb in bs.items():
    assert expand(X[i]*Y[j]-X[j]*Y[i]-baseL*baseB*A**3*(xx[i]-xx[j])*(p+A*bb))==0
for i in range(5):
    actual=X[i]**2+Y[i]**2
    target=prod(p+A*b for e,b in bs.items() if i in e)
    defect=Poly(expand(actual-target),p)
    assert defect.degree()<=1 and defect.nth(1)==A**3*E5*xx[i]*(E2*xx[i]-S*E2-E3)
    print('mark',i,'linear coefficient / A³',expand(defect.nth(1)/A**3))

# Rational parametrization of a dominant family of leading directions.
from sympy import cancel
b0,c0,k0,m0=symbols('b0 c0 k0 m0')
a0=b0*c0*(b0+c0)/(k0**2-b0*c0)
s0=a0+b0+c0
g0=a0*b0+a0*c0+b0*c0
u0=k0*a0
assert cancel(u0**2-s0*a0*b0*c0)==0
d0=(g0*s0-2*u0*m0)/(m0**2-g0)
# Identity before substitution, keeping this certificate small.
ss,gg,uu,mm=symbols('ss gg uu mm')
dd=(gg*ss-2*uu*mm)/(mm**2-gg)
assert cancel((uu+mm*dd)**2-gg*dd*(dd+ss)-uu**2)==0
