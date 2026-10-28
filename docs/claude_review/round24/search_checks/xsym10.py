import sympy as sp
m,Z=sp.symbols('m Z')
A=m**5-9*m**4+18*m**3+54*m**2+13*m-13
B=-4*m**5+36*m**4-104*m**3+72*m**2+108*m-108
D=-m**5+9*m**4-18*m**3-54*m**2+243*m-243
R=A*Z**4+B*(Z**2-Z)+D
print('R irreducible over Q:',sp.Poly(R,m,Z).is_irreducible)
disc=sp.factor(sp.discriminant(R,Z))
print('disc_Z R =',disc)
# Galois group of the generic quartic: resolvent cubic of A Z^4 + B Z^2 - B Z + D (monic after dividing by A)
p=B/A; q=-B/A; r=D/A
y=sp.symbols('y')
res=sp.factor(sp.together(y**3-p*y**2-4*r*y+(4*p*r-q**2)))
print('resolvent cubic factor:',res)
# specialisations: Galois group check at a few rational m
for mv in [sp.Rational(5,2),sp.Rational(7,3),2,sp.Rational(9,4)]:
    poly=sp.Poly(R.subs(m,mv),Z)
    print('m =',mv,' factorization over Q:',sp.factor_list(poly.as_expr()))
