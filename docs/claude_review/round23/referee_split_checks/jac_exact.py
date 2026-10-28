"""Exact Jacobian rank of (a,b,c,t) -> (F_j/F_0)_{j=1..4} on the conic z^2 = m x^2 + B, rationally parametrised
through (a, z_a). Upgrades split.md Thm 2(v) from numerical to exact."""
import sympy as sp, itertools
a,b,c,t = sp.symbols('a b c t')
d = a+b+c
sq=[a**2,b**2,c**2,d**2]; S=sum(sq)
e3=sum(p*q*r for p,q,r in itertools.combinations(sq,3)); m=a*b*c*d; B=(e3-S*m)/4
za=(a**2*b + a**2*c + a*b**2 + a*c**2 - b**2*c - b*c**2)/2
x=(a*(t**2+m)-2*za*t)/(t**2-m); z=za+t*(x-a)
assert sp.simplify(z**2-(m*x**2+B))==0
F0=(x-a)*(x-b)*(x-c)*(x+d); F1=x**4-S/2*x**2+m+2*x*z
F2=(x**2-b**2)*(x**2-c**2); F3=(x**2-a**2)*(x**2-c**2); F4=(x**2-a**2)*(x**2-b**2)
R=[F/F0 for F in (F1,F2,F3,F4)]
vars_=(a,b,c,t)
J=sp.Matrix([[sp.diff(r,v) for v in vars_] for r in R])
# homogeneity check: invariance under (a,b,c,t)->(l a, l b, l c, l^2 t)
l=sp.symbols('l')
for r in R:
    assert sp.simplify(r.subs({a:l*a,b:l*b,c:l*c,t:l**2*t},simultaneous=True)-r)==0
print("invariance under (a,b,c,t)->(la,lb,lc,l^2 t): OK  => generic rank <= 3")
for pt in [(1,2,4,sp.Rational(3,7)),(1,3,5,sp.Rational(-2,5)),(2,3,7,sp.Rational(11,3)),(sp.Rational(1,2),5,6,sp.Rational(1,9))]:
    Jn=J.subs(dict(zip(vars_,pt)))
    Jn=sp.Matrix(Jn.applyfunc(sp.nsimplify))
    print(pt,"exact rank:",Jn.rank())
