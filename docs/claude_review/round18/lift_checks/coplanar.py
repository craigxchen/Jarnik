"""Symbolic check (in x1..x4): aux_1, aux_2, aux_3 (S = empty over y_1,y_2,y_3) and B_bc (over y_4) are coplanar,
i.e. 3H - aux1 - aux2 - aux3 - B_bc ~ 2H = T^*O(1) identically on the family; plus the pencil identity."""
import sympy as sp
x1, x2, x3, x4, z, T, lam = sp.symbols('x1 x2 x3 x4 z T lam')
X = [x1, x2, x3, x4]
e1 = sum(X); e2 = sum(X[i]*X[j] for i in range(4) for j in range(i+1, 4))
e3 = sum(X[i]*X[j]*X[k] for i in range(4) for j in range(i+1,4) for k in range(j+1,4)); e4 = x1*x2*x3*x4
def PQ(l): return ((1 + l*e2)*z**2 - l*e3*z + l*e4, 1 + l*e1*z - l*z**2)
lb = 4*e4/(e3**2 - 4*e2*e4); lc = -4/e1**2
def A(l, Tv):
    P, Qq = PQ(l); poly = sp.Poly(sp.expand(P - Tv*Qq), z)
    c = poly.all_coeffs()[::-1] + [0]*3
    return c[0], c[1], c[2]
def point(j, sa, bpick, cpick):
    """sa=+1/-1 for a=+-x_j; bpick/cpick: 'x' for root x_j, 'o' for other root."""
    xj = X[j]; Tv = xj**2; a = sa*xj
    A0, A1, A2 = A(lb, Tv); b = xj if bpick == 'x' else sp.simplify(-A1/A2 - xj)
    u = sp.simplify((2*A2*b + A1)/a)
    C0, C1, C2 = A(lc, Tv); c = xj if cpick == 'x' else sp.simplify(-C1/C2 - xj)
    w = sp.simplify(2*C2*c + C1)
    return [a, u, w, 1]
rows = [point(0, -1, 'o', 'o'), point(1, -1, 'o', 'o'), point(2, -1, 'o', 'o'), point(3, -1, 'x', 'x')]
M = sp.Matrix(rows)
d = sp.factor(sp.together(M.det()))
print('det(aux1, aux2, aux3, B_bc) =', d)
# pencil identity P_lam - z^2 Q_lam = lam prod(z - x_j)
P, Qq = PQ(lam)
print('pencil identity:', sp.expand(P - z**2*Qq - lam*sp.prod([z - x for x in X])) == 0)
