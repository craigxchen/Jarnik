from sympy import *
from itertools import combinations
a,b,c,d,k,m,t=symbols('a b c d k m t')
ssum=a+b+c; g=a*b+a*c+b*c
y=[a,b,c,d,-ssum-d]
e4=sum(prod(cc) for cc in combinations(y,4))
assert expand(-e4-(g*d*(d+ssum)+ssum*a*b*c))==0
print("-e4(y) = g d(d+s) + s a b c : OK")
# identity sum 1/x = e4/e5: x=1/y, s4(x)=0 iff sum y=0 ; -s1 s5 = -e4(y)/e5(y)^2
X=symbols('X1:6'); Y=[1/xx for xx in X]
el=lambda v,r: sum(prod(cc) for cc in combinations(v,r))
assert simplify(el(X,4)/el(X,5)-sum(Y))==0
assert simplify(-el(X,1)*el(X,5)-(-el(Y,4)/el(Y,5)**2))==0
print("x=1/y dictionary OK")
# parametrization -> point in A^4 (a,b,c,d); Jacobian rank
A0=b*c*(b+c)/(k**2-b*c)
S0=A0+b+c; G0=A0*b+A0*c+b*c; u0=k*A0
D0=(G0*S0-2*u0*m)/(m**2-G0)
assert simplify(-(el([A0,b,c,D0,-S0-D0],4)) - (u0+m*D0)**2)==0
print("parametrized point satisfies -e4 = v^2 with v=u0+m d : OK")
F=Matrix([A0,b,c,D0]); Jac=F.jacobian([b,c,k,m])
vals={b:Rational(2),c:Rational(3),k:Rational(7),m:Rational(5,3)}
print("Jacobian rank at a random rational point:",Jac.subs(vals).rank())
# projective shape map x -> [x_i-x_5] from {(a,b,c,d)}: x=1/y, then differences; check image dimension 3
yv=[A0,b,c,D0,-S0-D0]; xv=[1/yy for yy in yv]
sh=[ (xv[i]-xv[4])/(xv[0]-xv[4]) for i in range(1,4)]
Js=Matrix(sh).jacobian([b,c,k,m])
print("rank of shape map (P^3 affine chart) Jacobian:",Js.subs(vals).rank())
# translation argument: e4(z+t) leading coefficient 5 in t
z=symbols('z1:6')
print("e4(z+t) leading coeff in t:",Poly(expand(el([zz+t for zz in z],4)),t).LC())
