"""Referee check (ref2_tri): symbolic verification of Astra II.5 identities."""
import sympy as sp

a, b, s, w, u, v, e = sp.symbols('a b s w u v e')  # e = sign epsilon (e^2=1)

# (II.1) from the determinant data: P1=(s,w), P2=(X2,Y2), |P1|^2=ab, |P2|^2=ac, det=+-a u
X2, Y2, dl = sp.symbols('X2 Y2 dl')
sub_ab = {s**2: a*b - w**2}
# det(P1,P2) = s*Y2 - w*X2 = dl*a*u  -> X2 = (s*Y2 - dl*a*u)/w, with Y2 = e*v
X2expr = (s*e*v - dl*a*u)/w
c_expr = sp.expand((X2expr**2 + v**2)/a)
c_expr = sp.expand(c_expr.subs(s**2, a*b - w**2)).subs({e**2: 1, dl**2: 1})
target = (a*u**2 + b*v**2 - 2*e*dl*s*u*v)/w**2
print('II.1 from determinant data:', sp.simplify(c_expr - target) == 0)

# Heron via the three-vector identity
A, B, C = a*u**2, b*v**2, None
c = (a*u**2 + b*v**2 + 2*e*s*u*v)/w**2
C = c*w**2
U = (u*v*w)**2
heron = sp.expand((A + B - C)**2 - (4*A*B - 4*U))
heron = sp.expand(heron.subs(e**2, 1).subs(s**2, a*b - w**2))
print('Heron (A+B-C)^2=4AB-4U :', sp.simplify(heron) == 0)
print('AB-U = s^2u^2v^2       :', sp.simplify(sp.expand(A*B - U - s**2*u**2*v**2).subs(s**2, a*b-w**2) - 0) == 0
      or sp.simplify(sp.expand(A*B-U).subs(a*b, s**2+w**2)))

# reverse parametrization identities
r1 = sp.expand(w**2*(a*c - v**2) - (a*u + e*s*v)**2)
r1 = sp.expand(sp.expand(r1).subs(e**2, 1).subs(s**2, a*b - w**2))
r2 = sp.expand(w**2*(b*c - u**2) - (s*u + e*b*v)**2)
r2 = sp.expand(sp.expand(r2).subs(e**2, 1).subs(s**2, a*b - w**2))
print('ac-v^2 = ((au+esv)/w)^2:', sp.simplify(r1) == 0)
print('bc-u^2 = ((su+ebv)/w)^2:', sp.simplify(r2) == 0)

# vector identity det(Y,Z)X+det(Z,X)Y+det(X,Y)Z=0
x1, x2, y1, y2, z1, z2 = sp.symbols('x1 x2 y1 y2 z1 z2')
det = lambda p, q: p[0]*q[1] - p[1]*q[0]
X, Y, Z = (x1, x2), (y1, y2), (z1, z2)
vec = [sp.expand(det(Y, Z)*X[k] + det(Z, X)*Y[k] + det(X, Y)*Z[k]) for k in range(2)]
print('det(Y,Z)X+det(Z,X)Y+det(X,Y)Z=0:', vec == [0, 0])
aa = [det(Y, Z)*X[k] for k in range(2)]
bb = [det(Z, X)*Y[k] for k in range(2)]
print('det of first two = prod of dets:', sp.expand(det(aa, bb) - det(X, Y)*det(Y, Z)*det(Z, X)) == 0)

# Vieta: other root of C^2-2(A+B)C+(A-B)^2+4U
Cs, As, Bs, Us = sp.symbols('C A B U')
poly = sp.expand((As + Bs - Cs)**2 - 4*As*Bs + 4*Us)
roots = sp.solve(poly, Cs)
print('Vieta roots:', roots, ' product:', sp.simplify(roots[0]*roots[1]))
