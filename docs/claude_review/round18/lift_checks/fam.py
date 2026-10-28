"""Rational 3-parameter family of balanced genus-one curves in Mbar_{0,6}.

Parameters: distinct nonzero rationals x1..x4 (common arguments); anchors x1,x2,x3.
f_a(z) = z^2.  Pencil: f_lam = P_lam/Q_lam,
    P_lam = (1 + lam e2) z^2 - lam e3 z + lam e4,   Q_lam = 1 + lam e1 z - lam z^2,
so that P_lam(x_j) = x_j^2 Q_lam(x_j) for j=1..4 (since P - z^2 Q = lam * prod(z - x_j)).
lam_b = 4 e4/(e3^2 - 4 e2 e4)  : other member with branch value 0;
lam_c = -4/e1^2                : other member with branch value infinity.
Exact arithmetic only (fractions.Fraction)."""
from fractions import Fraction as Fr

def esym(xs):
    e = [Fr(1)]
    for x in xs:
        e = [e[0]] + [e[k] + x * e[k - 1] for k in range(1, len(e))] + [x * e[-1]]
    return e  # e[0]=1, e[1]=sum, ..., coefficients of prod (z + x)... careful sign below

def family(xs):
    xs = [Fr(x) for x in xs]
    e = esym(xs)  # elementary symmetric e1..e4 (prod (1 + x t) coefficients)
    e1, e2, e3, e4 = e[1], e[2], e[3], e[4]
    def PQ(lam):
        P = [lam * e4, -lam * e3, 1 + lam * e2]   # low -> high
        Q = [Fr(1), lam * e1, -lam]
        return P, Q
    lam_b = 4 * e4 / (e3 ** 2 - 4 * e2 * e4)
    lam_c = -4 / e1 ** 2
    return dict(xs=xs, e=(e1, e2, e3, e4), maps={'a': PQ(Fr(0)), 'b': PQ(lam_b), 'c': PQ(lam_c)},
                lam={'a': Fr(0), 'b': lam_b, 'c': lam_c})

def ev(p, z):
    r = Fr(0)
    for c in reversed(p):
        r = r * z + c
    return r

def dev(p, z):
    return ev([k * p[k] for k in range(1, len(p))], z)

def fval(PQ, z):
    P, Q = PQ
    return ev(P, z) / ev(Q, z)

def fder(PQ, z):
    P, Q = PQ
    q = ev(Q, z)
    return (dev(P, z) * q - ev(P, z) * dev(Q, z)) / q ** 2

def other_root(PQ, T, z0):
    """Other root of P(z) - T Q(z) = 0 given root z0."""
    P, Q = PQ
    A = [P[k] - T * Q[k] for k in range(3)]  # A0 + A1 z + A2 z^2
    if A[2] == 0:
        return None  # root at infinity
    return -A[1] / A[2] - z0

def disc_T(PQ):
    """Discriminant of P(z) - T Q(z) as polynomial in T (low->high)."""
    P, Q = PQ
    # (A1)^2 - 4 A2 A0 with Ak = P[k] - T Q[k]
    def lin(k):
        return [P[k], -Q[k]]
    def mulp(u, v):
        r = [Fr(0)] * (len(u) + len(v) - 1)
        for i, a in enumerate(u):
            for j, b in enumerate(v):
                r[i + j] += a * b
        return r
    A1sq = mulp(lin(1), lin(1)); A20 = mulp(lin(2), lin(0))
    return [A1sq[i] - 4 * A20[i] for i in range(3)]

def check(xs, verbose=True):
    F = family(xs)
    xs = F['xs']; maps = F['maps']
    ok = True
    # common arguments
    for j, x in enumerate(xs):
        for k in 'abc':
            P, Q = maps[k]
            if ev(Q, x) == 0 or fval(maps[k], x) != x * x:
                ok = False; print('common-argument failure', k, x)
    # branch values
    D = {k: disc_T(maps[k]) for k in 'abc'}
    if verbose:
        print('disc_T a,b,c:', [[str(c) for c in D[k]] for k in 'abc'])
    # transversality: derivatives at each common argument pairwise distinct, finite, nonzero
    for j, x in enumerate(xs):
        ds = [fder(maps[k], x) for k in 'abc']
        if len(set(ds)) < 3 or any(d == 0 for d in ds):
            ok = False; print('transversality failure at x_%d' % (j + 1), ds)
    return F, D, ok

if __name__ == '__main__':
    import sys
    xs = [Fr(s) for s in sys.argv[1:5]] if len(sys.argv) >= 5 else [Fr(2), Fr(3), Fr(-5), Fr(1)]
    F, D, ok = check(xs)
    print('lam_b, lam_c =', F['lam']['b'], F['lam']['c'], ' ok =', ok)
