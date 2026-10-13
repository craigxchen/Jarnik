from fractions import Fraction as Fr
import qpoly as Q
from ratrec import nullspace
def reconstruct(f, D):
    pts = []
    k = 0
    while len(pts) < 2 * D + 8:
        x = Fr(3 * k + 1, 11) + 2; k += 1
        try:
            pts.append((x, Fr(f(x))))
        except ZeroDivisionError:
            pass
    rows = [[x ** i for i in range(D + 1)] + [-y * x ** j for j in range(D + 1)] for x, y in pts]
    ns = nullspace(rows)
    v = ns[0]
    p = Q.norm(v[:D + 1]); q = Q.norm(v[D + 1:])
    g = Q.gcd(p, q)
    p, _ = Q.divmod_(p, g); q, _ = Q.divmod_(q, g)
    lc = q[-1]
    p, q = Q.scal(1 / lc, p), Q.scal(1 / lc, q)
    nchk = 0
    for k in range(40):
        x = Fr(13 * k + 5, 17) - 7
        try:
            fx = f(x)
        except ZeroDivisionError:
            continue
        assert Q.ev(p, x) == fx * Q.ev(q, x), 'reconstruction check failed'
        nchk += 1
        if nchk >= 6: break
    return p, q
