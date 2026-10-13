# Referee: a single relation on a rank-3 lattice has NONDEGENERATE zeros arbitrarily close to the identity.
# mu = integer measure on Z^3 with ord(mu) = 2 (kills 1, y1, y2, y3).  ghat(s) = sum mu_y e^{i<y,s>}
#   = -Q(s)/2 - i C(s)/6 + O(|s|^4),  Q(s) = sum mu_y <y,s>^2,  C(s) = sum mu_y <y,s>^3.
# If Q and C have a common zero u on S^2 with transversal intersection, then G(t,v) = (Re ghat(tv)/t^2,
# Im ghat(tv)/t^3) is real-analytic near (0,u) (Re ghat is even in s, Im ghat odd), G(0,u)=0 and D_vG(0,u) is
# invertible on T_u S^2, so the IFT gives zeros s(t) = t v(t) of ghat, v(0)=u, for all small t.  If u lies in no
# rational plane, then for each nonzero integer y, t -> <y,s(t)> is analytic and not identically 0, so all but
# countably many t give NONDEGENERATE zeros.  All algebra below is exact (sympy).
import sympy as sp
pts = [(0,0,0),(1,0,0),(0,1,0),(0,0,1),(1,1,0),(2,0,1),(1,2,-1),(0,-1,2)]
Mm = sp.Matrix([[1,a,b,c] for (a,b,c) in pts]).T
ns = Mm.nullspace()
# integer combination
v = ns[0]
for k in range(1, len(ns)): v = v + (k+1)*ns[k]
den = sp.ilcm(*[sp.fraction(e)[1] for e in v]); v = v*den
mu = {p: int(v[i]) for i, p in enumerate(pts)}
print("mu =", mu)
s1, s2, s3, X, Y = sp.symbols('s1 s2 s3 X Y')
lin = lambda p: p[0]*s1 + p[1]*s2 + p[2]*s3
for d in range(2): assert sp.expand(sum(c*lin(p)**d for p, c in mu.items())) == 0
Q = sp.expand(sum(c*lin(p)**2 for p, c in mu.items()))
C = sp.expand(sum(c*lin(p)**3 for p, c in mu.items()))
print("Q =", Q); print("C =", C)
# affine chart s1 = 1
q = Q.subs({s1: 1, s2: X, s3: Y}); cc = C.subs({s1: 1, s2: X, s3: Y})
res = sp.Poly(sp.resultant(q, cc, Y), X)
print("resultant in X:", res.as_expr())
fac = sp.factor_list(res.as_expr())
print("factorisation:", fac)
for f, e in fac[1]:
    pf = sp.Poly(f, X)
    if pf.degree() < 1: continue
    for rt in sp.real_roots(pf):
        xv = sp.N(rt, 30)
        # find Y common root
        ys = [yv for yv in sp.Poly(q.subs(X, xv), Y).nroots(n=30) if abs(sp.im(yv)) < 1e-20]
        for yv in ys:
            yv = sp.re(yv)
            if abs(sp.N(cc.subs({X: xv, Y: yv}))) < 1e-15:
                J = sp.Matrix([[sp.diff(q, X), sp.diff(q, Y)], [sp.diff(cc, X), sp.diff(cc, Y)]]).subs({X: xv, Y: yv})
                print(f"common real zero direction u=(1, {sp.N(xv,15)}, {sp.N(yv,15)}), factor degree {pf.degree()}, "
                      f"multiplicity {e}, Jacobian det = {sp.N(J.det(),10)}")
# ---- exact certification that u is irrational (no rational plane) and transversal ----
mp = sp.Poly(5488*X**4 + 5504*X**3 - 4032*X**2 + 256*X + 497, X)
assert mp.is_irreducible
# remainder of cc by q as polynomials in Y with coefficients in Q[X]
r1 = sp.prem(sp.Poly(cc, Y), sp.Poly(q, Y))          # pseudo-remainder, degree <= 1 in Y
r1 = sp.Poly(r1, Y)
A_, B_ = r1.coeff_monomial(Y), r1.coeff_monomial(1)
K = sp.QQ.algebraic_field  # not used; reduce via polynomial arithmetic mod mp
def red(e):
    return sp.rem(sp.Poly(sp.expand(e), X), mp)
Ared, Bred = red(A_), red(B_)
# inverse of Ared mod mp
inv = sp.invert(Ared.as_expr(), mp.as_expr(), X)
yexpr = red(-Bred.as_expr()*inv)
print("y = polynomial in x (mod minimal poly):", yexpr.as_expr())
coeffs = sp.Poly(yexpr.as_expr(), X).all_coeffs()[::-1] + [0]*4
print("coefficients of x^2, x^3 in y:", coeffs[2], coeffs[3], "-> (1,x,y) Q-linearly independent:",
      (coeffs[2] != 0 or coeffs[3] != 0))
# check q, cc vanish exactly at (x, y(x)) modulo mp
print("q(x,y(x)) mod mp =", red(q.subs(Y, yexpr.as_expr())).as_expr(), "; cc(x,y(x)) mod mp =", red(cc.subs(Y, yexpr.as_expr())).as_expr())
Jdet = sp.Matrix([[sp.diff(q, X), sp.diff(q, Y)], [sp.diff(cc, X), sp.diff(cc, Y)]]).det()
Jr = red(Jdet.subs(Y, yexpr.as_expr()))
print("Jacobian det mod mp (nonzero element of Q(x) iff transversal):", Jr.as_expr())
# numeric illustration: follow the zero curve s(t)=t v(t) by Newton for a few t, show |ghat| ~ 0 and no small relation
import mpmath as mpm
mpm.mp.dps = 60
def ghat(s):
    return sum(c*mpm.e**(1j*(p[0]*s[0]+p[1]*s[1]+p[2]*s[2])) for p, c in mu.items())
xv = [r for r in sp.Poly(mp.as_expr(), X).nroots(n=50) if abs(sp.im(r)) < 1e-40]
for x0 in xv:
    x0 = mpm.mpf(str(sp.re(x0)))
    y0 = mpm.mpf(str(sp.N(yexpr.as_expr().subs(X, sp.Float(str(x0), 60)), 60)))
    u = mpm.matrix([1, x0, y0]); u = u/mpm.norm(u)
    for t in [mpm.mpf('1e-2'), mpm.mpf('1e-4'), mpm.mpf('1e-8')]:
        # solve ghat(t*(u + a*e1 + b*e2)) = 0 for (a,b) with e1,e2 basis of u-perp
        e1 = mpm.matrix([-u[1], u[0], 0]); e1 = e1/mpm.norm(e1)
        e2 = mpm.matrix([u[1]*e1[2]-u[2]*e1[1], u[2]*e1[0]-u[0]*e1[2], u[0]*e1[1]-u[1]*e1[0]])
        def F(a, b):
            s = t*(u + a*e1 + b*e2)
            g = ghat(s)
            return [mpm.re(g)/t**2, mpm.im(g)/t**3]
        a, b = mpm.findroot(F, (mpm.mpf(0), mpm.mpf(0)))
        s = t*(u + a*e1 + b*e2)
        print(f"  t={mpm.nstr(t,3)}: |ghat(s)| = {mpm.nstr(abs(ghat(s)),3)}, |s| = {mpm.nstr(mpm.norm(s),3)}, "
              f"s/|s| = {[mpm.nstr(z,8) for z in (s/mpm.norm(s))]}")
