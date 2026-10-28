import sympy as sp, pickle
b1,b2,X,Y,T2=sp.symbols('b1 b2 X Y T2')
b0=sp.Integer(1); b3=-1-b1-b2
bs=[b0,b1,b2,b3]
T2v=sp.expand(b0*b1*b2*b3)
# n = prod (X^2 - 2 i b X - (t^2/b^2 + b^2)); split real/imag with X real: A_r = X^2 - |p|^2, B_r = -2 b X
A=[X**2-(T2v/b**2+b**2) for b in bs]; B=[-2*b*X for b in bs]
# product of (A_r + i B_r)
re,im=sp.Integer(1),sp.Integer(0)
for a_,b_ in zip(A,B):
    re,im=sp.expand(re*a_-im*b_),sp.expand(re*b_+im*a_)
sig=-sp.prod([ (T2v/b**2+b**2) for b in bs])
L=sp.together(re-im**2/(4*sig)+sig)
num,den=sp.fraction(L)
num=sp.expand(num)
PX=sp.Poly(num,X)
assert all(PX.coeff_monomial(X**k)==0 for k in range(1,9,2))
ell=sp.Poly(sum(PX.coeff_monomial(X**(2*k))*Y**k for k in range(5)),Y)
print('ell(0) numerator =',sp.simplify(ell.coeff_monomial(1)))
cub=sp.Poly(sp.expand(sp.cancel(ell.as_expr()/Y)),Y)
print('cubic degree',cub.degree())
disc=sp.discriminant(cub.as_expr(),Y)
disc=sp.factor(disc)
pickle.dump((disc,cub.as_expr(),den),open('xsym4.pkl','wb'))
fl=sp.factor_list(disc)
print('number of factors',len(fl[1]))
for f,e in fl[1]:
    print('  factor deg',sp.Poly(f,b1,b2).total_degree(),'mult',e, ' nterms',len(sp.Poly(f,b1,b2).terms()))
