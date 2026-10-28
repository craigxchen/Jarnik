"""General Pell-balanced family in split LM_5 (generalising m-czech (a,b,c)=(1,2,4)).
roots r in {a,b,c,d}, d=a+b+c; n=prod(x^2-r^2); S=sum r^2, e2,e3 elementary symmetric of the squares;
m=abcd; A=(S^2/4+2m-e2)/4, B=(e3-S m)/4; z^2=A x^2+B; f5=x^4-(S/2)x^2+m+2xz.
Checks (exact, symbolic in a,b,c): N(f5)=n on the conic; A r^2+B is a square at each root; profile (B.C, cuts) at
random rational (a,b,c); orders at infinity; Jacobian rank of (a,b,c,x) -> torus point (mod scaling) at random points."""
import sympy as sp, random
from fractions import Fraction as Fr
a, b, c, x, z = sp.symbols('a b c x z')
d = a + b + c
R = [a, b, c, d]
sq = [r**2 for r in R]
S = sum(sq)
e2 = sum(sq[i]*sq[j] for i in range(4) for j in range(i+1, 4))
e3 = sum(sq[i]*sq[j]*sq[k] for i in range(4) for j in range(i+1, 4) for k in range(j+1, 4))
m = a*b*c*d
A = sp.expand((S**2/4 + 2*m - e2)/4); B = sp.expand((e3 - S*m)/4)
n = sp.expand(sp.prod([x**2 - s for s in sq]))
f5 = x**4 - S/2*x**2 + m + 2*x*z
f5b = x**4 - S/2*x**2 + m - 2*x*z
N5 = sp.expand(sp.expand(f5*f5b).subs(z**2, A*x**2 + B))
print("N(f5) - n == 0 identically:", sp.simplify(N5 - n) == 0)
print("A =", sp.factor(A)); print("B =", sp.factor(B))
for r in R:
    val = sp.factor(A*r**2 + B)
    print("  A r^2 + B at r =", r, ":", val)
D0 = (x - a)*(x - b)*(x - c)*(x + d)
P = lambda p, q: (x**2 - p**2)*(x**2 - q**2)
F = [D0, f5, P(b, c), P(a, c), P(a, b)]
for j in range(1, 5):
    print("  F%d - F0 =" % j, sp.collect(sp.expand(F[j] - F[0]), x))
# profile at random rational parameters
random.seed(1)
def profile(av, bv, cv):
    sub = {a: av, b: bv, c: cv}
    Av, Bv = A.subs(sub), B.subs(sub)
    places = []
    for r in [av, bv, cv, av+bv+cv]:
        for s in (r, -r):
            w = sp.sqrt(Av*s**2 + Bv)
            for zz in (w, -w): places.append((s, zz))
    vecs = []
    for (s, zz) in places:
        vecs.append(tuple(1 if sp.simplify(f.subs(sub).subs({x: s, z: zz})) == 0 else 0 for f in F))
    BC = sum(max(v) - min(v) for v in vecs)
    cuts = set()
    for v in vecs:
        if max(v) != min(v):
            Aset = frozenset(j for j in range(5) if v[j]); Bset = frozenset(range(5)) - Aset
            cuts.add(min(tuple(sorted(Aset)), tuple(sorted(Bset))))
    return BC, len(cuts), Av, Bv, len(set(places))
for (av, bv, cv) in [(1, 2, 4), (1, 3, 5), (2, 3, 7), (1, 5, 6), (Fr(1), Fr(3, 2), Fr(5, 3))]:
    av, bv, cv = sp.Rational(av), sp.Rational(bv), sp.Rational(cv)
    print("(a,b,c)=", (av, bv, cv), "-> B.C, #cuts, A, B, #places:", profile(av, bv, cv))
# Jacobian rank of (a,b,c,x) -> (log F_j - log F_0)_{j=1..4}, z = +sqrt(A x^2+B)
zz = sp.sqrt(A*x**2 + B)
G = [f.subs(z, zz)/F[0] for f in F[1:]]
J = sp.Matrix([[sp.diff(g, v) for v in (a, b, c, x)] for g in G])
for t in range(4):
    pt = {a: sp.Rational(random.randint(1, 9), random.randint(1, 5)), b: sp.Rational(random.randint(1, 9), 3),
          c: sp.Rational(random.randint(10, 30), 7), x: sp.Rational(random.randint(40, 90), 1)}
    import numpy as np
    Jn = np.array(J.subs(pt).evalf(30).tolist(), dtype=float)
    sv = sorted(np.linalg.svd(Jn, compute_uv=False))
    print("Jacobian singular values at", pt, ":", ["%.3e" % s for s in sv])
