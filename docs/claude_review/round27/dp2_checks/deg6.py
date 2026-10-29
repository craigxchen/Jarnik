"""Curves of class -3K (degree 6) through P with a 4-fold point at P.
Sections of -3K: cubic forms in x,y,z (10) and s*linear (3); local expansion at P with s = w(x,y) (z=1).
Order >= 4 at P: 10 linear conditions -> 3-dim space (a net), containing x*D and y*D.
alpha < 2 for an irreducible curve C in |-3K| through P needs (MR Thm 2.16 + D.C >= 0) a 4-fold point whose
tangent cone is kappa*q^2 (two conjugate double branches, alpha = 6/4) or l^4 (one rational 4-fold branch, 6/4).
We find all members of the net whose tangent cone is a square up to scalar, and test rationality."""
import sympy as sp
from sympy import symbols, expand, Poly, Matrix, Rational
from geom import default_f, local_branch
x, y, z, s = symbols('x y z s')
f = default_f()
w = local_branch(f, order=6)
wser = sum(w)  # Taylor polynomial of w up to degree 5

mons = [x**i*y**j*z**(3-i-j) for i in range(4) for j in range(4-i)] + [s*x, s*y, s*z]
def local(expr):
    e = expand(expr.subs({z: 1}).subs({s: wser}))
    P = Poly(e, x, y)
    return {mon: co for mon, co in P.terms() if sum(mon) <= 5}
locs = [local(mo) for mo in mons]
low = [(i, j) for i in range(4) for j in range(4 - i)]  # monomials of degree <= 3 in x,y
M = Matrix([[loc.get(mon, 0) for loc in locs] for mon in low])
ker = M.nullspace()
print('dim of sections of -3K with order >= 4 at P:', len(ker))
secs = [expand(sum(v[i]*mons[i] for i in range(len(mons)))) for v in ker]
def tc(sec):
    loc = local(sec)
    return expand(sum(co*x**i*y**j for (i, j), co in loc.items() if i + j == 4))
TCs = [tc(S_) for S_ in secs]
for T in TCs:
    print('tangent cone:', sp.factor(T))
# general member: l0*sec0 + l1*sec1 + l2*sec2 ; tangent cone = sum l_i TC_i ; require = (q0 x^2 + q1 xy + q2 y^2)^2
l0, l1, l2, q0, q1, q2 = symbols('l0 l1 l2 q0 q1 q2')
TC = expand(l0*TCs[0] + l1*TCs[1] + l2*TCs[2])
Q2 = expand((q0*x**2 + q1*x*y + q2*y**2)**2)
eqs = Poly(expand(TC - Q2), x, y).coeffs()
# normalise one of l's to 1 (try each chart)
sols_all = []
for chart in [l0, l1, l2]:
    eqc = [e.subs(chart, 1) for e in eqs]
    unknowns = [v for v in (l0, l1, l2, q0, q1, q2) if v != chart]
    G = sp.groebner(eqc, *unknowns, order='lex')
    print('chart', chart, 'Groebner basis length', len(G.exprs), 'last element degree',
          Poly(G.exprs[-1], unknowns[-1]).degree() if G.exprs[-1].free_symbols <= {unknowns[-1]} else 'n/a')
    last = G.exprs[-1]
    if last.free_symbols <= {unknowns[-1]}:
        print('  eliminant factorisation:', sp.factor_list(last))
    sols_all.append((chart, G))
