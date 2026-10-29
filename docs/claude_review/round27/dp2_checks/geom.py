"""Exact geometry of a degree-2 del Pezzo surface S: w^2 = f(x,y,z) at P = (0:0:1:1).
Usage: python3 geom.py   (prints certificate data for the fixed choice of f below)
All computations exact (sympy rationals, flint integer polynomials)."""
import sys, itertools
import sympy as sp
from sympy import Poly, symbols, expand, Rational, factor_list, resultant, discriminant, gcd
import flint

x, y, z, t, m, c = symbols('x y z t m c')

def default_f():
    # the surface S of round 27 lane dp2 (search_f.py seed 1); f(0,0,1) = 1
    return (x**3*y + 2*x**2*y*z - 2*x*y**3 - 2*x*y*z**2 + 2*x*z**3 + 2*y**4 + 2*y**3*z
            + y**2*z**2 - y*z**3 + z**4)

def old_singular_f():
    # first trial (singular at (-1:0:1)); kept for the record
    return (z**4 + x*z**3 + (y**2 - x*y)*z**2 + (x**3 - x**2*y + 2*y**3)*z
            + (x**4 - x**2*y**2 + x*y**3 + 3*y**4))

def parts(f):
    """homogeneous parts g_i (degree i in x,y) of f(x,y,1)"""
    F = Poly(expand(f.subs(z, 1)), x, y)
    g = [0]*5
    for (i, j), co in F.terms():
        g[i+j] += co*x**i*y**j
    return [expand(gi) for gi in g]

def local_branch(f, order=5):
    """Taylor polynomial of w = sqrt(f(x,y,1)) at (0,0), w(0,0)=+1, by degree"""
    g = parts(f)
    assert g[0] == 1
    eps = symbols('eps')
    G = sum(g[i]*eps**i for i in range(1, 5))
    ser = sp.series(sp.sqrt(1 + G), eps, 0, order).removeO()
    ser = expand(ser)
    return [expand(ser.coeff(eps, i)) for i in range(order)]

def analyse(f, verbose=True):
    out = {}
    w = local_branch(f)
    T2hat = expand(z**2 + z*w[1] + w[2])     # homogenised 2-jet: section of -2K
    F4 = expand(f - T2hat**2)                 # pi(D) = {F4 = 0}
    P4 = Poly(F4, x, y, z)
    # check triple point at p=(0:0:1): all monomials have degree >= 3 in (x,y)
    mindeg = min(i + j for (i, j, k), co in P4.terms())
    out['F4_min_xy_degree'] = mindeg
    F3 = expand(sum(co*x**i*y**j for (i, j, k), co in P4.terms() if i + j == 3))
    F4p = expand(sum(co*x**i*y**j for (i, j, k), co in P4.terms() if i + j == 4))
    out['F3'] = F3; out['F4prime'] = F4p
    out['w3'] = w[3]
    out['F3_equals_2w3'] = expand(F3 - 2*w[3]) == 0
    # binary cubic F3 in slope: F3(1,t) and check root at infinity
    F3t = Poly(expand(F3.subs({x: 1, y: t})), t)
    out['F3_lead_y3'] = Poly(F3, x, y).coeff_monomial(y**3)
    out['F3_factor'] = factor_list(F3)
    out['disc_F3'] = discriminant(F3t)
    out['F3_real_roots'] = len(sp.Poly(F3t).real_roots())
    # D irreducible iff Res(F3, F4') != 0 (no line through p is a component)
    F4pt = Poly(expand(F4p.subs({x: 1, y: t})), t)
    out['Res_F3_F4p'] = resultant(F3t, F4pt)
    # also check common root at infinity (x=0): coefficients of y^3 in F3, y^4 in F4'
    out['inf_check'] = (Poly(F3, x, y).coeff_monomial(y**3), Poly(F4p, x, y).coeff_monomial(y**4))
    return out, w, T2hat, F4

if __name__ == '__main__':
    f = default_f()
    out, w, T2hat, F4 = analyse(f)
    print('f =', f)
    print('w1 =', w[1]); print('w2 =', w[2]); print('w3 =', w[3])
    print('T2hat =', T2hat)
    for k, v in out.items():
        print(k, ':', v)
