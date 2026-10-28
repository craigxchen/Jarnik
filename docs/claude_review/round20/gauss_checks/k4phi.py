import sympy as sp
u1,v1,u2,v2,u3,v3 = sp.symbols('u1 v1 u2 v2 u3 v3', real=True)
r = [u1+sp.I*v1, u2+sp.I*v2, u3+sp.I*v3]
s1 = sum(r); s2 = r[0]*r[1]+r[0]*r[2]+r[1]*r[2]; s3 = r[0]*r[1]*r[2]
I = lambda e: sp.im(sp.expand(e)); R_ = lambda e: sp.re(sp.expand(e))
Phi = sp.expand(I(s3)*I(s1) + (R_(s1)*I(s1) - I(s2))*I(s2) - I(s1)**2*R_(s2))
print("Phi =", sp.factor(Phi))
# specialise rho = i (u1=0,v1=1)
P = sp.factor(Phi.subs({u1:0, v1:1}))
print("Phi(i,a,c) =", P)
print("degree in (u3,v3):", sp.Poly(Phi.subs({u1:0,v1:1}), u3, v3).total_degree())
