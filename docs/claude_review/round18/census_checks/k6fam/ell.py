"""Quotient elliptic curve E': W^2 = Delta_b(T) Delta_c(T) of the fiber product, exact group law,
sections over the four common values, non-torsion via Mazur (order of a torsion point in {1..10,12})."""
from fractions import Fraction as Fr
from math import isqrt
from fam import family, disc_T, check

def polymul(u, v):
    r = [Fr(0)] * (len(u) + len(v) - 1)
    for i, a in enumerate(u):
        for j, b in enumerate(v):
            r[i + j] += a * b
    while len(r) > 1 and r[-1] == 0:
        r.pop()
    return r

def ev(p, z):
    r = Fr(0)
    for c in reversed(p):
        r = r * z + c
    return r

def is_sq(q):
    q = Fr(q)
    if q < 0: return None
    a, b = q.numerator, q.denominator
    ra, rb = isqrt(a), isqrt(b)
    return Fr(ra, rb) if ra * ra == a and rb * rb == b else None

class Curve:
    """y^2 = x^3 + A2 x^2 + A4 x + A6 (general a1=a3=0 form)."""
    def __init__(self, A2, A4, A6):
        self.A2, self.A4, self.A6 = Fr(A2), Fr(A4), Fr(A6)
    def on(self, P):
        if P is None: return True
        x, y = P
        return y * y == x ** 3 + self.A2 * x * x + self.A4 * x + self.A6
    def add(self, P, Q):
        if P is None: return Q
        if Q is None: return P
        x1, y1 = P; x2, y2 = Q
        if x1 == x2 and y1 == -y2:
            return None
        if P == Q:
            m = (3 * x1 * x1 + 2 * self.A2 * x1 + self.A4) / (2 * y1)
        else:
            m = (y2 - y1) / (x2 - x1)
        x3 = m * m - self.A2 - x1 - x2
        y3 = m * (x1 - x3) - y1
        return (x3, y3)
    def mul(self, n, P):
        R = None
        for _ in range(n):
            R = self.add(R, P)
        return R

def quotient_curve(xs):
    F = family(xs)
    Db = disc_T(F['maps']['b']); Dc = disc_T(F['maps']['c'])
    g = polymul(Db, Dc)  # cubic in T (Dc linear)
    assert len(g) == 4 and g[0] == 0, g
    k = g[3]
    # W^2 = k T^3 + g2 T^2 + g1 T ; X = k T, Y = k W  ->  Y^2 = X^3 + g2 X^2 + k g1 X
    C = Curve(g[2], k * g[1], 0)
    pts = []
    for x in F['xs']:
        T = x * x
        w = is_sq(ev(g, T))
        pts.append(None if w is None else (k * T, k * w))
    return F, C, g, pts

def nontorsion(C, P):
    if P is None: return False
    for n in list(range(1, 11)) + [12]:
        if C.mul(n, P) is None:
            return False
    return True

if __name__ == '__main__':
    import sys
    xs = [Fr(s) for s in sys.argv[1:5]] if len(sys.argv) >= 5 else [Fr(2), Fr(3), Fr(-5), Fr(1)]
    F, D, ok = check(xs, verbose=False)
    F, C, g, pts = quotient_curve(xs)
    print('params', [str(x) for x in xs], 'transversal/nondeg ok:', ok)
    print("E': Y^2 = X^3 + (%s) X^2 + (%s) X" % (C.A2, C.A4))
    for j, P in enumerate(pts):
        print('  section over x_%d^2:' % (j + 1), None if P is None else (str(P[0]), str(P[1])), ' on curve:', C.on(P), ' non-torsion (Mazur):', nontorsion(C, P))
