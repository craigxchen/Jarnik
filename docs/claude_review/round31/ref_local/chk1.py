import sympy as sp
m,Y,W=sp.symbols('m Y W')
c=m**3-11*m**2+39*m-13; D=m**4+70*m**2-48*m+9
u=-(m+3)*c/(16*(m-3)**2*(m-1)**2); v=-(m+1)*c/(4*(m-3)**3*(m-1))
G=sp.expand((Y**2-Y+u)**2+v*Y)
print("v/u =",sp.factor(v/u))
print("u(m)-u(-m) =",sp.factor(u-u.subs(m,-m)))
# discriminant of monic quartic via generic formula
a3,a2,a1,a0=sp.symbols('a3 a2 a1 a0')
gen=Y**4+a3*Y**3+a2*Y**2+a1*Y+a0
dg=sp.discriminant(gen,Y)
P=sp.Poly(G,Y); co=P.all_coeffs(); assert co[0]==1
dG=dg.subs({a3:co[1],a2:co[2],a1:co[3],a0:co[4]})
claim=u*(2*D)*(2*(m+1)*c/((m-3)**5*(m-1)**3))**2
print("disc identity:",sp.simplify(sp.together(dG-claim))==0)
print("factor disc:",sp.factor(dG))
# b-quartic and Y relation: roots b of W^4+vW^2-v^2W+uv^2; Y=-b^2/v satisfies Y^2-Y+u=b?
bq=W**4+v*W**2-v**2*W+u*v**2
Yb=-W**2/v
r=sp.rem(sp.numer(sp.together(Yb**2-Yb+u-W)),sp.numer(sp.together(bq)),W)
print("Y^2-Y+u-b mod bquartic:",sp.simplify(r))
print("res check G(Y(b)) mod bq:",sp.simplify(sp.rem(sp.numer(sp.together(G.subs(Y,Yb))),sp.numer(sp.together(bq)),W)))
print("limits u,v at oo:",sp.limit(u,m,sp.oo),sp.limit(v,m,sp.oo))
print("u(0)=",u.subs(m,0))
print("D irreducible:",sp.factor_list(D), "c:",sp.factor_list(c))
# m^2 recovery: solve v/u=k
k=sp.symbols('k'); print(sp.solve(sp.Eq(4*(m**2-1),k*(m**2-9)),m))
