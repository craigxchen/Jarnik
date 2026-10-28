import sympy as sp
t,b1,b2,b3,X=sp.symbols('t b1 b2 b3 X',real=True)
b0=-(b1+b2+b3)
bs=[b0,b1,b2,b3]; p=[t/b+sp.I*b for b in bs]
mb=lambda z: -sp.conjugate(z)
roots={1:[p[1],p[0],mb(p[2]),mb(p[3])],2:[p[2],p[0],mb(p[1]),mb(p[3])],3:[p[3],p[0],mb(p[1]),mb(p[2])]}
for j in (1,2,3):
    P=sp.expand(sp.prod([X-r for r in roots[j]]))
    ims=[sp.factor(sp.together(sp.im(sp.expand(P.coeff(X,k))))) for k in (3,2,1)]
    print('j',j,ims)
