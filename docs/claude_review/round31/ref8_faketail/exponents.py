# Exact-rational re-derivation of the tail exponents of VIII.8 and VIII.10.
from fractions import Fraction as Fr
def F(g,nu,d,b,atom):   # per-token base-q exponent: count + atoms + X_k + gap
    return (1-g)*(1-nu) + (b - atom(g))*nu - d/2
# VIII.8
d=Fr(199,100); b=Fr(1708,1000); a=Fr(12,100)
strong=lambda g: a+Fr(1,2)+g/2; weak=lambda g: a+Fr(1,2)
# F is linear in nu and in g on each piece; check vertices and the nu-coefficient sign
for name,atom,glo,ghi in (("strong",strong,Fr(8,100),Fr(10,100)),("weak",weak,Fr(10,100),Fr(1))):
    vals=[F(g,nu,d,b,atom) for g in (glo,ghi) for nu in (Fr(0),Fr(1,2))]
    print("VIII.8",name,"max over box =",max(vals),"=",float(max(vals)))
print(" closed forms: strong 49/1000-3g/4 at g=.08:",Fr(49,1000)-Fr(3,4)*Fr(8,100),
      " weak 49/1000-g/2 at g=.1:",Fr(49,1000)-Fr(1,2)*Fr(1,10))
print(" nonprimitive: 1/2+b/4-d/2 =",Fr(1,2)+b/4-d/2,"; zeta: b/4*(q)-d/2*(q/2) per q =",b/4-d/4)
print(" B/(HE): b-.08-1.625 =",b-Fr(8,100)-Fr(1625,1000),"; P/B vs H^2:",d-b,"vs",Fr(16,100))
print(" E vs qA sqrt q: 1.625 vs",1+a+Fr(1,2))
# VIII.10: symbolic check of the exponent per token at nu=1/2
import sympy as sp
xi,ka,ta,aa,ep,g,nu=sp.symbols('xi kappa tau a epsilon gamma nu')
dd=2-xi; bb=sp.Rational(3,2)+aa+ka+2*ep
Fs=sp.expand((1-g)*(1-nu)+(bb-(aa+sp.Rational(1,2)+g/2))*nu-dd/2)
Fw=sp.expand((1-g)*(1-nu)+(bb-(aa+sp.Rational(1,2)))*nu-dd/2)
print("VIII.10 dF/dnu strong:",sp.diff(Fs,nu)," weak:",sp.diff(Fw,nu))
print(" F at nu=1/2 strong:",sp.simplify(Fs.subs(nu,sp.Rational(1,2)))," weak:",sp.simplify(Fw.subs(nu,sp.Rational(1,2))))
print(" strong at g=kappa:",sp.simplify(Fs.subs({nu:sp.Rational(1,2),g:ka})),"; weak at g=tau:",sp.simplify(Fw.subs({nu:sp.Rational(1,2),g:ta})))
print(" nonprimitive (1-kappa)/2+b/4-d/2 with b->d:",sp.simplify((1-ka)/2+dd/4-dd/2))
# closure: lower endpoints kappa=2xi, tau=3xi, a=3xi in 1.5+a+kappa<2-xi
print(" closure:",sp.solve(sp.Eq(sp.Rational(3,2)+3*xi+2*xi,2-xi),xi))
# explicit feasible point for xi=1/12 - 1/1000
x=Fr(1,12)-Fr(1,1000); k=2*x+Fr(1,10000); t=k+x+Fr(1,10000); A_=t+Fr(1,10000)
e=Fr(1,10**6); B_=Fr(3,2)+A_+k+2*e
print(" feasible at xi=",x,":",B_<2-x, k>2*x+4*e, t>k+x+2*e, A_>t, t<Fr(1,2), "margin d-b=",float(2-x-B_))
print(" limits: d->",2-Fr(1,12),"=",Fr(23,12)," 4d->",4*Fr(23,12))
print(" VIII.9: 1.01/d at d=23/12:",float(Fr(101,100)/Fr(23,12)),"; middle case 4/3-.53/2=",float(Fr(4,3)-Fr(53,200)),
      "; s>=q/2 case (16/9-.53)/2+4/9=",float((Fr(16,9)-Fr(53,100))/2+Fr(4,9)),"; second 1-.53/2=",1-Fr(53,200))
