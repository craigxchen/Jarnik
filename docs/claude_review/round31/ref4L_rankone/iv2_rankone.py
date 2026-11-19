# IV.2 family: delta_i as a polynomial in x_i reduced mod prod(X-x_j) (coeffs in s1..s5,P,q).
# Decide when delta_i = c * x_i^2 (the only possible rank-one shape, l(V_i) = L*B*x_i).
import sympy as sp
X,P,q,s1,s2,s3,s5,l=sp.symbols('X P q s1 s2 s3 s5 l')
s4=0
f=X**5-s1*X**4+s2*X**3-s3*X**2+s4*X-s5      # prod (X-x_j)
# b_ij = x_i x_j (s1-x_i-x_j); n_i=prod_{j!=i}(P+q b_ij): symmetric fn of the other roots.
# Use: prod_{j!=i}(T - x_j) = f(T)/(T-x_i); b_ij as polynomial in x_j: y -> x y (s1-x-y) = -x y^2 + x(s1-x) y
# n_i = prod_{j!=i} (P + q*(-x y_j^2 + x(s1-x) y_j)) = resultant-type; compute via resultant in y.
y=sp.Symbol('y')
g=sp.quo(f.subs(X,y)-f.subs(X,X), y-X)   # prod_{j!=i}(y-x_j) with x_i=X  (monic quartic in y)
g=sp.expand(sp.cancel((f.subs(X,y)-f)/(y-X)))
h=P+q*(-X*y**2+X*(s1-X)*y)
# prod over roots y_j of g of h(y_j) = Res_y(g,h) (g monic)
n=sp.resultant(g,h,y)
n=sp.expand(n)
b0=sp.Rational(3,2)*s5
B=b0*q**2; L=l*q
hX=X**3-s1*X**2+s2*X
U=P+q*hX
V1=P*U+B*X; V2=L*U
delta=sp.expand(V1**2+V2**2-n)
delta=delta.subs(l**2,-s1*s5)
r=sp.rem(sp.Poly(delta,X),sp.Poly(f,X))
r=sp.Poly(sp.expand(r.as_expr()),X)
for k in range(5):
    print('X^%d:'%k, sp.factor(r.coeff_monomial(X**k)))
