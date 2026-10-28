"""The pencil lifts C'_lam = {f_a(a) = f_b(b) = f_c(c) = f_lam(d)} of E_t: exact checks at t, lam.
Class beta_7: at each y_j the 16 points (+-x_j, {x_j,b'}, {x_j,c'}, {x_j,d'}) give 15 (j<=3) resp. 11 (j=4)
contacts, transversal iff the slopes f_m'(x_j) (m = a,b,c,lam) are nonzero and pairwise distinct, and the other
roots differ from x_j.  Genus: (Z/2)^4-cover of P^1_T branched at {0,inf,s,t'} u {branch values of f_lam};
genus 9 iff the two branch values of f_lam are new (not in {0, inf, s, t'}) and distinct; 5 if exactly one is new."""
import sys
from fractions import Fraction as Fr
from fam import family, ev, dev, fder, disc_T, check

def pencil(F, lam):
    e1, e2, e3, e4 = F['e']
    return ([lam * e4, -lam * e3, 1 + lam * e2], [Fr(1), lam * e1, -lam])

def run(xs, lam):
    F, D, ok = check(xs, verbose=False); assert ok
    xs = F['xs']; lam = Fr(lam)
    maps = dict(F['maps']); maps['d'] = pencil(F, lam)
    # branch values of each map: roots of disc_T
    s = -D['b'][1] / D['b'][2]; tp = -D['c'][0] / D['c'][1]
    Dd = disc_T(maps['d'])
    okslope = True
    for x in xs:
        sl = [fder(maps[m], x) for m in 'abcd']
        if any(v == 0 for v in sl) or len(set(sl)) < 4:
            okslope = False
    # branch values of f_lam: roots of Dd (quadratic in T)
    import sympy as sp
    Tt = sp.symbols('T')
    poly = sp.Poly(sum(sp.Rational(str(c)) * Tt ** i for i, c in enumerate(Dd)), Tt)
    roots = sp.roots(poly, multiple=True) if poly.degree() > 0 else []
    old = [0, sp.Rational(str(s)), sp.Rational(str(tp))]
    new = [r for r in roots if all(sp.simplify(r - o) != 0 for o in old)]
    deg = poly.degree()
    genus = None
    nb = len(set(new)) + (0 if deg == 2 else 1)   # degree drop = branch value at infinity (old)
    if deg == 2 and len(set(roots)) == 2:
        genus = {2: 9, 1: 5, 0: 1}[len(new)]
    print('t =', [str(x) for x in xs], ' lam =', lam, ' slopes ok:', okslope, ' disc_T(f_lam) degree', deg,
          ' branch values', [str(r) for r in roots], ' new:', len(new), ' genus:', genus)
    return okslope, genus

if __name__ == '__main__':
    for t in ([2, 3, -5, 1], [Fr(1, 3), 3, 2, Fr(5, 3)]):
        for lam in (1, 2, -1, Fr(1, 7)):
            run([Fr(x) for x in t], lam)
