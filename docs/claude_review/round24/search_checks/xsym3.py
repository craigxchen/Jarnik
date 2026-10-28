import sympy as sp
b1,b2,X,Y,T2=sp.symbols('b1 b2 X Y T2',real=True)
b0=sp.Integer(1); b3=-1-b1-b2
bs=[b0,b1,b2,b3]
T2v=b0*b1*b2*b3   # t^2
n=sp.expand(sp.prod([X**2-2*sp.I*b*X-(T2/b**2+b**2) for b in bs]))
n=sp.expand(n.subs(T2,T2v))
Ren=sp.expand(sp.re(n)); Imn=sp.expand(sp.im(n))
P_Im=sp.Poly(Imn,X); P_Re=sp.Poly(Ren,X)
print('Im n degree',P_Im.degree(),' coeffs:',[sp.factor(c) for c in P_Im.all_coeffs()][:3])
c5=sp.factor(sp.together(P_Im.coeff_monomial(X**5)))
print('[X^5] Im n =',c5)
sig=-sp.prod([ (T2v/b**2+b**2) for b in bs])
L=sp.expand(Ren-Imn**2/(4*sig)+sig)
print('L(0) =',sp.simplify(L.subs(X,0)))
