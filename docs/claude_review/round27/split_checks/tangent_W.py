# Tangent directions of the W-curves at e (round 23 Thm 2(iii)): (Y1..Y4) ~ (+-2 sqrt(abcd), (a+b)(a+c), (a+b)(b+c), (a+c)(b+c)),
# d = a+b+c. Implicit equation of the surface of directions Sigma_W in P^3 and its irreducibility.
import sympy as sp
a,b,c,Y1,Y2,Y3,Y4,l=sp.symbols('a b c Y1 Y2 Y3 Y4 l')
d=a+b+c
p,q,r=a+b,a+c,b+c
G=4*Y1**2*Y2**2*Y3**2*Y4**2-(2*(Y2**2+Y3**2+Y4**2)*Y2**2*Y3**2*Y4**2-(Y2**4*Y3**4+Y2**4*Y4**4+Y3**4*Y4**4))
sub={Y1:2*sp.sqrt(a*b*c*d),Y2:p*q,Y3:p*r,Y4:q*r}
print("vanishes on directions:", sp.simplify(G.subs(sub))==0)
print("factorisation:", sp.factor_list(sp.expand(G)))
