# Judge check J1 (independent code): Theorem B ingredients and form-independence.
import sympy as sp
m,W,Y,lam=sp.symbols('m W Y lam')
c = m**3-11*m**2+39*m-13
u = -(m+3)*c/(16*(m-3)**2*(m-1)**2)
v = -(m+1)*c/(4*(m-3)**3*(m-1))
bq = W**4 + v*W**2 - v**2*W + u*v**2
# (1) Y = -W^2/v satisfies G(Y) = (Y^2-Y+u)^2 + v Y : resultant check
G = (Y**2-Y+u)**2 + v*Y
res = sp.resultant(sp.numer(sp.together(bq)), sp.numer(sp.together(Y + W**2/v)), W)
q = sp.cancel(sp.together(res)/sp.together(G))
print("Res_W(bquartic, Y+W^2/v)/G is free of Y:", not sp.together(q).has(Y))
# (2) scale invariants: roots b with e1=0,e2=v,e3=v^2,e4=u v^2 ; u=e4/e2^2, v=e3^2/e2^3, invariant under b->lam b
e2,e3,e4=sp.symbols('e2 e3 e4')
U = e4/e2**2; V = e3**2/e2**3
print("u invariant under b->lam b:", sp.simplify(U.subs({e2:lam**2*e2,e3:lam**3*e3,e4:lam**4*e4},simultaneous=True)-U)==0)
print("v invariant under b->lam b:", sp.simplify(V.subs({e2:lam**2*e2,e3:lam**3*e3,e4:lam**4*e4},simultaneous=True)-V)==0)
print("check u,v from b-quartic coefficients:", sp.simplify(U.subs({e2:v,e3:v**2,e4:u*v**2})-u)==0, sp.simplify(V.subs({e2:v,e3:v**2,e4:u*v**2})-v)==0)
# (3) discriminant of G in Y, factor, odd-order finite points
Gn = sp.Poly(sp.numer(sp.together(G)), Y)
D = sp.factor(sp.discriminant(Gn))
print("disc_Y(numerator of G) =", D)
num = sp.numer(sp.together(D)); 
fl = sp.factor_list(num)
odd = [(f,e) for f,e in fl[1] if e%2==1 and f.has(m)]
print("odd-exponent factors:", odd, " #odd points =", sum(sp.degree(f,m) for f,e in odd))
lc = sp.Poly(sp.numer(sp.together(G)),Y).LC()
print("leading coeff in Y:", sp.factor(lc))
for f,e in odd:
    print("  gcd with lc:", sp.gcd(f, lc))
# distinct odd points: squarefree product
P = sp.Mul(*[f for f,e in odd])
print("odd-point polynomial squarefree:", sp.degree(sp.gcd(P, sp.diff(P,m)),m)==0)
# behaviour at infinity: degree of disc numerator vs denominator parity, done via m=1/s
s=sp.symbols('s')
Dinf = sp.factor(sp.together(D.subs(m,1/s)))
ns = sp.numer(Dinf); ds = sp.denom(Dinf)
print("order at m=inf (s=0): numer", sp.Poly(ns,s).terms()[-1][0], "denom", sp.Poly(ds,s).terms()[-1][0] if ds.has(s) else 0)
# (4) Galois group of G_m over Q at sample m
for mv in [sp.Rational(-2), sp.Rational(-3,2), sp.Rational(5), sp.Rational(7,3), sp.Rational(-5,4)]:
    g = sp.Poly(sp.numer(sp.together(G.subs(m,mv))), Y)
    try:
        grp = sp.polys.numberfields.galoisgroups.galois_group(g, by_name=True)
        print("m=",mv," irreducible:", g.is_irreducible, " Gal:", grp)
    except Exception as ex:
        print("m=",mv,"galois_group failed:", ex, "irreducible:", g.is_irreducible)
