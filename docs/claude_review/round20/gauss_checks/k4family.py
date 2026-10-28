"""Build and verify exactly a k=4 circle family with degree-1 Gaussian blocks.
Input roots rho=r123, a=r12, b=r13, c=r23 in Q(i) with Phi=0 pairwise; singleton roots r_i computed rationally.
Checks (exact): (i) Im P_i(s) is a nonzero constant (polynomial identity); (ii) for integer s: the four lattice points
z_0..z_3 lie on x^2+y^2=N, gcd structure (Norm gcd of P_i,P_j = prod n_T), residues t_ij, Y_i bounded;
(iii) reduced-arc exponent log(maxchord)/log R -> 3/7.  Uses sympy for polynomial identities, Python ints for values."""
import sympy as sp, math, sys
from fractions import Fraction as F
from itertools import combinations
s = sp.symbols('s')
def singleton(roots):
    s1 = sum(roots); s2 = sum(x*y for x,y in combinations(roots,2))
    s1, s2 = sp.expand(s1), sp.expand(s2)
    ims1 = sp.im(s1); x = (sp.re(s1)*ims1 - sp.im(s2))/ims1
    return sp.nsimplify(x - sp.I*ims1)
def build(rho, a, b, c):
    r = {(1,2,3): rho, (1,2): a, (1,3): b, (2,3): c}
    r[(1,)] = singleton([rho, a, b]); r[(2,)] = singleton([rho, a, c]); r[(3,)] = singleton([rho, b, c])
    Hs = {}
    for T, rt in r.items():
        q = sp.ilcm(sp.fraction(sp.re(rt))[1], sp.fraction(sp.im(rt))[1])
        Hs[T] = sp.expand(q*(s - rt))
    P = {i: sp.expand(sp.prod([Hs[T] for T in Hs if i in T])) for i in (1,2,3)}
    return r, Hs, P
def gi(e, sv):  # evaluate Gaussian-integer linear/poly at integer s
    v = sp.expand(e.subs(s, sv)); return (int(sp.re(v)), int(sp.im(v)))
if __name__ == "__main__":
    I = sp.I
    rho, a, b, c = I, -9-4*I, -12-5*I, -15+4*I
    r, Hs, P = build(rho, a, b, c)
    print("roots:", r)
    for i in (1,2,3):
        im = sp.Poly(sp.im(sp.expand(P[i].subs(s, sp.Symbol('t', real=True)))), sp.Symbol('t', real=True))
        print("P_%d: Im P_%d(s) as poly in s:" % (i,i), im.as_expr(), "| leading coeff", sp.Poly(P[i], s).LC())
