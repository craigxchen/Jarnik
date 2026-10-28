import sympy as sp
a0,a1,a2,a3,b1,b2,b3,X=sp.symbols('a0 a1 a2 a3 b1 b2 b3 X',real=True)
b0=-(b1+b2+b3)
p=[a0+sp.I*b0,a1+sp.I*b1,a2+sp.I*b2,a3+sp.I*b3]
mb=lambda z: -sp.conjugate(z)
roots={1:[p[1],p[0],mb(p[2]),mb(p[3])],2:[p[2],p[0],mb(p[1]),mb(p[3])],3:[p[3],p[0],mb(p[1]),mb(p[2])]}
eqs=[]
for j in (1,2,3):
    P=sp.expand(sp.prod([X-r for r in roots[j]]))
    for k in (3,2,1):
        c=sp.expand(P.coeff(X,k)); eqs.append(sp.simplify(sp.im(c)))
for e in eqs: print(sp.factor(e))
