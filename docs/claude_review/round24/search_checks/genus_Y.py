# genus of C1': (Y^2 - Y + u(m))^2 + v(m) Y = 0, Y = b^2/(-e2) (scale- and sign-invariant), degree-4 cover of the m-line
import sympy as sp
m,Y=sp.symbols('m Y')
cub=m**3-11*m**2+39*m-13
u=-(m+3)*cub/(16*(m-3)**2*(m-1)**2); v=-(m+1)*cub/(4*(m-3)**3*(m-1))
G=sp.expand((Y**2-Y+u)**2+v*Y)
num,den=sp.fraction(sp.together(G))
Gp=sp.Poly(sp.expand(num),Y)
print('leading coeff in Y:',sp.factor(Gp.LC()))
d=sp.factor(sp.discriminant(sp.expand(num),Y))
print('disc_Y =',d)
print('irreducible over Q(m)?',sp.Poly(sp.expand(num),Y,m).is_irreducible)
fl=sp.factor_list(sp.expand(num),Y)
for f,e in fl[1]:
    print('factor (deg in Y = %d):'%sp.Poly(f,Y).degree(), sp.factor(f))
