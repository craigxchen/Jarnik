import sympy as sp, pickle
u,v,k,e2,e3,e4=sp.symbols('u v k e2 e3 e4')
F=pickle.load(open('F24.pkl','rb'))
Phi=sp.expand(F.subs({e4:u*e2**2,e3:sp.sqrt(v*e2**3)})/e2**12)
Phi=sp.expand(sp.simplify(Phi))
print('Phi =',Phi)
G=sp.expand(Phi.subs(v,k*u)/u**3)
P=sp.Poly(G,u); print('degree in u:',P.degree())
c2,c1,c0=[P.coeff_monomial(u**d) for d in (2,1,0)]
D=sp.factor(sp.expand(c1**2-4*c2*c0))
print('discriminant in k:',D)
print('c0 (cubic part):',sp.factor(c0)); print('c2:',sp.factor(c2))
