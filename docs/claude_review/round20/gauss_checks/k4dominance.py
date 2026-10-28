"""Exact Jacobian ranks (forward-mode AD over Q, exact Fractions) for the k=4 degree-1 circle family.
Variables (u2,v2,u3,v3,u4,v4,s) = (Re,Im of r12, r13, r23, s); r123 = i fixed (normalisation by real affine maps);
constraints Phi(i,a,b)=Phi(i,a,c)=Phi(i,b,c)=0 (3 eqs).  Singleton roots are rational functions (sing()).
Targets: (1) similarity class: tan arg(z_x/z_0), x=1,2,3 (coordinates on the real points of T_4 = norm-one torus mod
diagonal, birational to M^tw_{0,6}); (2) P(X_4) mod real scaling: tan arg z_x, x=0..3.
Solution set is a smooth 4-manifold near the point iff rank d(constraints) = 3; the map to T_4 is a submersion iff
rank [dcons; dtargets1] = 6, to P(X_4) iff rank [dcons; dtargets2] = 7.
Also: the rational 2-parameter subfamily u2=u3=u4=u, x_j=1+v_j, x_a+x_b+x_c=1, u^2=x_a x_b x_c."""
from fractions import Fraction as F
import sympy as sp
NV = 7
class D:
    def __init__(s, v, g=None): s.v = F(v); s.g = g if g is not None else [F(0)]*NV
    def __add__(a, b): b = b if isinstance(b, D) else D(b); return D(a.v+b.v, [x+y for x,y in zip(a.g,b.g)])
    __radd__ = __add__
    def __neg__(a): return D(-a.v, [-x for x in a.g])
    def __sub__(a, b): return a + (-(b if isinstance(b, D) else D(b)))
    def __rsub__(a, b): return D(b) - a
    def __mul__(a, b): b = b if isinstance(b, D) else D(b); return D(a.v*b.v, [a.v*y+b.v*x for x,y in zip(a.g,b.g)])
    __rmul__ = __mul__
    def __truediv__(a, b):
        b = b if isinstance(b, D) else D(b); q = a.v/b.v
        return D(q, [(x - q*y)/b.v for x,y in zip(a.g,b.g)])
class C:
    def __init__(s, re, im): s.re = re if isinstance(re, D) else D(re); s.im = im if isinstance(im, D) else D(im)
    def __add__(a, b): return C(a.re+b.re, a.im+b.im)
    def __sub__(a, b): return C(a.re-b.re, a.im-b.im)
    def __mul__(a, b): return C(a.re*b.re - a.im*b.im, a.re*b.im + a.im*b.re)
    def conj(a): return C(a.re, -a.im)
def var(i, val): g = [F(0)]*NV; g[i] = F(1); return D(val, g)
def phi(a, b):
    (x2,y2),(x3,y3) = a, b
    return y2*y3*(x2-x3)*(x2-x3) + y2*x2*x2 + y3*x3*x3 + (1+y2)*(1+y3)*(y2+y3)*(1+y2+y3)
def sing(t):
    rs = [C(x, y) for x, y in t]
    s1 = rs[0]+rs[1]+rs[2]; s2 = rs[0]*rs[1] + rs[0]*rs[2] + rs[1]*rs[2]
    x = (s1.re*s1.im - s2.im)/s1.im
    return C(x, -s1.im)
def evaluate(vals):
    u2,v2,u3,v3,u4,v4,s = [var(i, vals[i]) for i in range(NV)]
    A,B,Cc = (u2,v2),(u3,v3),(u4,v4); rho = (D(0), D(1))
    cons = [phi(A,B), phi(A,Cc), phi(B,Cc)]
    roots = {(1,2,3): C(0,1), (1,2): C(u2,v2), (1,3): C(u3,v3), (2,3): C(u4,v4),
             (1,): sing([rho,A,B]), (2,): sing([rho,A,Cc]), (3,): sing([rho,B,Cc])}
    H = {T: C(s, D(0)) - r for T, r in roots.items()}
    z = []
    for x in range(4):
        acc = C(1, 0)
        for T in H: acc = acc * (H[T] if x in T else H[T].conj())
        z.append(acc)
    rel = [(z[x]*z[0].conj()) for x in (1,2,3)]
    return cons, [w.im/w.re for w in rel], [w.im/w.re for w in z]
def rank(rows): return sp.Matrix([[sp.Rational(x.numerator, x.denominator) for x in r] for r in rows]).rank()
if __name__ == "__main__":
    for pt in [(-6,-4,-6,-3,-6,5,7), (-6,-4,-6,-3,-6,5,F(13,2)), (-8,-3,-10,-6,F(-27,2),8,5),
               (F(-3,2),F(-5,2),F(-3,2),F(-3,2),F(-3,2),2,11)]:
        cons, t1, t2 = evaluate(pt)
        assert all(c.v == 0 for c in cons), [c.v for c in cons]
        Jc = [c.g for c in cons]
        print("point", [str(x) for x in pt], ": rank dcons =", rank(Jc), "; rank [dcons; T_4] =", rank(Jc + [t.g for t in t1]),
              "(6 => submersion onto T_4) ; rank [dcons; P(X_4)] =", rank(Jc + [t.g for t in t2]), "(7 => onto P(X_4))")
    def sub_pt(X, m, s):
        xb = X*(1-X)/(m*m+X); xc = 1 - X - xb; u = m*xb
        return (u, X-1, u, xb-1, u, xc-1, s)
    for (X, m, s) in [(F(-3), F(-3), 7), (F(-2), F(1,3), 9), (F(5,2), F(2), 4)]:
        pt = sub_pt(X, m, s); cons, t1, t2 = evaluate(pt)
        print("subfamily X=%s m=%s s=%s -> a,b,c =" % (X,m,s), [str(x) for x in pt[:6]],
              "constraints hold:", all(c.v == 0 for c in cons))
