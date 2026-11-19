"""Symbolic check of the (P2)->(P3) algebra and the perfect-matching identity in II.7."""
import sympy as sp, itertools
u = sp.symbols('u1:6'); lam = sp.symbols('lam')
c = {}
for i in range(5):
    for j in range(5):
        if i < j: c[(i,j)] = sp.Symbol('c%d%d'%(i+1,j+1))
for (i,j) in list(c): c[(j,i)] = -c[(i,j)]
t = {(i,j): lam*u[i]*u[j]*c[(i,j)] for (i,j) in c}
D = [sp.prod([c[(i,j)] for j in range(5) if j != i]) for i in range(5)]
U = sp.prod(u)
ok1 = all(sp.simplify(sp.prod([t[(i,j)] for j in range(5) if j!=i]) - lam**4*u[i]**3*U*D[i]) == 0 for i in range(5))
print("prod_j t_ij = lam^4 u_i^3 (prod u) D_i :", ok1)
ok2 = True
for i in range(5):
    o = [k for k in range(5) if k != i]
    for (a,b),(cc,d) in [((o[0],o[1]),(o[2],o[3])),((o[0],o[2]),(o[1],o[3])),((o[0],o[3]),(o[1],o[2]))]:
        m = (t[(a,b)]/c[(a,b)])*(t[(cc,d)]/c[(cc,d)])
        ok2 &= sp.simplify(m - lam**2*U/u[i]) == 0
print("m_i = lam^2 prod(u)/u_i for all three matchings:", ok2)
# minors of W_i=(u_i, u_i l_i) vs t_ij L_ij / lam, with l_i=<A_i,(p,q)> and c_ij L_ij=<A_j-A_i,(p,q)>
p,q = sp.symbols('p q'); A = [(sp.Symbol('a%d'%i), sp.Symbol('b%d'%i)) for i in range(5)]
l = [a*p+b*q for a,b in A]
ok3 = all(sp.expand(u[i]*u[j]*l[j]-u[j]*u[i]*l[i] - u[i]*u[j]*((A[j][0]-A[i][0])*p+(A[j][1]-A[i][1])*q))==0 for i in range(5) for j in range(5))
print("det(W_i,W_j) = u_i u_j <A_j-A_i,(p,q)> :", ok3)
