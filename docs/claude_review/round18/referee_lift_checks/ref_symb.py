import sympy as sp
x1,x2,x3,x4,lam,T,z=sp.symbols('x1 x2 x3 x4 lam T z')
X=[x1,x2,x3,x4]
e1=sum(X); e2=sum(X[i]*X[j] for i in range(4) for j in range(i+1,4))
e3=x1*x2*x3+x1*x2*x4+x1*x3*x4+x2*x3*x4; e4=x1*x2*x3*x4
lb=4*e4/(e3**2-4*e2*e4); lc=-4/e1**2
def A(l,TT): return (1+l*e2+l*TT, -l*(e3+e1*TT), l*e4-TT)
def Delta(l,TT):
    a2,a1,a0=A(l,TT); return a1**2-4*a2*a0
Db=sp.Poly(sp.together(Delta(lb,T)).as_numer_denom()[0],T); Dbd=sp.together(Delta(lb,T)).as_numer_denom()[1]
kb=sp.cancel(Db.coeff_monomial(T**2)/Dbd); s=sp.cancel(-Db.coeff_monomial(T)/Db.coeff_monomial(T**2))
assert sp.simplify(Db.coeff_monomial(1))==0
Dc=sp.together(Delta(lc,T)); Dcn,Dcd=Dc.as_numer_denom(); Dcn=sp.Poly(sp.expand(Dcn),T)
assert Dcn.degree()==1
kc=sp.cancel(Dcn.coeff_monomial(T)/Dcd); tp=sp.cancel(-Dcn.coeff_monomial(1)/Dcn.coeff_monomial(T))
def pt(j,a_sign,b_is_x,c_is_x):
    xj=X[j]; y=xj**2
    a2b,a1b,a0b=A(lb,y); a2c,a1c,a0c=A(lc,y)
    bp=sp.cancel(-a1b/a2b-xj); cp=sp.cancel(-a1c/a2c-xj)
    a=a_sign*xj; b=xj if b_is_x else bp; c=xj if c_is_x else cp
    u=sp.cancel((2*a2b*b+a1b)/a); w=sp.cancel(2*a2c*c+a1c)
    return [a,u,w,sp.Integer(1)]
aux=[pt(j,-1,False,False) for j in range(3)]
Bbc=pt(3,-1,True,True)
M=sp.Matrix(aux+[Bbc])
# exact symbolic determinant via rational function arithmetic
det=sp.cancel(sp.together(M.det(method='berkowitz')))
print('det(aux1,aux2,aux3,B_bc) symbolic =',det)
# check model equations symbolically at B_bc
for P in aux+[Bbc]:
    a,u,w,h=P
    print('on model:',sp.cancel(u**2-kb*(a**2-s)), sp.cancel(w**2-kc*(a**2-tp)))
