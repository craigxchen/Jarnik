"""Exact Jacobian of (s1..s4,t) -> angles arg(z_j/z_0), j=1..4, at rational points."""
import sys
from sympy import symbols, I, Matrix, expand, diff, Rational, im, re, simplify, nsimplify
t = symbols('t', real=True)
S = symbols('s1:5', real=True)
EDGES = [(0, 1), (0, 2), (1, 2), (0, 3)]
def conjp(e, subs):  # conjugate of polynomial with real symbols: replace I by -I
    return expand(e.subs(I, -I))
def marks(n):
    A = Matrix([[t**2+1, t], [t, 1]])**n
    G = []
    for s in S:
        w = A*Matrix([1, s]); G.append(expand(w[0]+I*w[1]))
    Gc = [expand(g.subs(I, -I)) for g in G]
    z = [(t+2*I)*Gc[0]*Gc[1]*Gc[2]*Gc[3]]
    for (j, k) in EDGES:
        e = (t-2*I)*G[j]*G[k]
        for l in range(4):
            if l not in (j, k): e *= Gc[l]
        z.append(e)
    return z
def jac(n, pt):
    z = marks(n)
    vars_ = list(S)+[t]
    sub = dict(zip(vars_, pt))
    zv = [complex_rat(zz.subs(sub)) for zz in z]
    rows = []
    for j in range(1, 5):
        row = []
        for v in vars_:
            a = expand(diff(z[j], v)).subs(sub); b = expand(diff(z[0], v)).subs(sub)
            q = expand(a/z[j].subs(sub) - b/z[0].subs(sub))
            row.append(im(q))
        rows.append(row)
    M = Matrix(rows)
    return M.rank(), M[:, :4].rank(), M[:, :4].det()
def complex_rat(x): return expand(x)
if __name__ == '__main__':
    for n in (1, 2):
        for pt in [(1, 2, 3, 5, 3), (Rational(1, 2), Rational(1, 3), Rational(2, 7), Rational(-5, 3), 7), (Rational(3, 4), Rational(-1, 5), Rational(7, 2), 0, 10)]:
            r5, r4, d = jac(n, pt)
            print('n=%d point(s1..s4,t)=%s rank(4x5)=%d rank(seeds 4x4)=%d det=%s' % (n, pt, r5, r4, d.evalf(6)))
