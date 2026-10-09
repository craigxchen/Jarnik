import sys, pickle
sys.path.insert(0, '..')
from fractions import Fraction as Fr
from polylib import *
from slice2 import degs

def coeff_polys_in_x(P):
    dx, dy = degs(P)
    return [norm([P.get((i, j), Fr(0)) for j in range(dy + 1)]) for i in range(dx + 1)]

def disc_cubic(c):  # c = [d, c, b, a] polys in y (a x^3 + b x^2 + c x + d)
    d, cc, b, a = c
    t = lambda *ps: prod(list(ps))
    D = add(add(add(scal(18, t(a, b, cc, d)), scal(-4, t(b, b, b, d))), add(t(b, b, cc, cc), scal(-4, t(a, cc, cc, cc)))), scal(-27, t(a, a, d, d)))
    return D

def sqfree_info(D):
    g = pgcd(D, deriv(D))
    return deg(D), deg(g)

if __name__ == '__main__':
    base, a, b, P = pickle.load(open(sys.argv[1], 'rb'))
    cx = coeff_polys_in_x(P)
    assert len(cx) == 4
    D = disc_cubic(cx)
    print('disc_x degree in y:', deg(D), ' gcd(D,D\') degree:', sqfree_info(D)[1])
    Py = {(j, i): c for (i, j), c in P.items()}
    cy = coeff_polys_in_x(Py)
    D2 = disc_cubic(cy)
    print('disc_y degree in x:', deg(D2), ' gcd degree:', sqfree_info(D2)[1])
