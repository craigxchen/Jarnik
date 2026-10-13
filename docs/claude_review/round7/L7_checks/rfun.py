"""Exact rational functions in one variable x over Q: RF(num, den) with qpoly coefficient lists,
and polynomials in a second variable s with RF coefficients."""
from fractions import Fraction as Fr
import qpoly as Q


class RF:
    __slots__ = ('n', 'd')

    def __init__(self, n, d=None):
        if isinstance(n, RF):
            self.n, self.d = n.n, n.d
            return
        if not isinstance(n, list):
            n = [Fr(n)]
        if d is None:
            d = [Fr(1)]
        n, d = Q.norm(n), Q.norm(d)
        if not d:
            raise ZeroDivisionError
        if not n:
            self.n, self.d = [], [Fr(1)]
            return
        g = Q.gcd(n, d)
        if Q.deg(g) > 0:
            n, _ = Q.divmod_(n, g)
            d, _ = Q.divmod_(d, g)
        lc = d[-1]
        self.n, self.d = Q.scal(1 / lc, n), Q.scal(1 / lc, d)

    @staticmethod
    def x():
        return RF([Fr(0), Fr(1)])

    def __add__(s, o):
        o = o if isinstance(o, RF) else RF(o)
        return RF(Q.add(Q.mul(s.n, o.d), Q.mul(o.n, s.d)), Q.mul(s.d, o.d))
    __radd__ = __add__

    def __neg__(s):
        return RF(Q.scal(-1, s.n), s.d)

    def __sub__(s, o):
        return s + (-(o if isinstance(o, RF) else RF(o)))

    def __rsub__(s, o):
        return RF(o) - s

    def __mul__(s, o):
        o = o if isinstance(o, RF) else RF(o)
        return RF(Q.mul(s.n, o.n), Q.mul(s.d, o.d))
    __rmul__ = __mul__

    def __truediv__(s, o):
        o = o if isinstance(o, RF) else RF(o)
        if not o.n:
            raise ZeroDivisionError
        return RF(Q.mul(s.n, o.d), Q.mul(s.d, o.n))

    def __rtruediv__(s, o):
        return RF(o) / s

    def iszero(s):
        return not s.n

    def ev(s, x):
        return Q.ev(s.n, x) / Q.ev(s.d, x)

    def __repr__(s):
        return 'RF(deg %d / deg %d)' % (Q.deg(s.n), Q.deg(s.d))


def pnorm(p):
    p = list(p)
    while p and p[-1].iszero():
        p.pop()
    return p


def padd(p, q):
    n = max(len(p), len(q))
    z = RF(0)
    return pnorm([(p[i] if i < len(p) else z) + (q[i] if i < len(q) else z) for i in range(n)])


def psub(p, q):
    return padd(p, [-c for c in q])


def pscal(c, p):
    return pnorm([c * a for a in p])


def pmul(p, q):
    if not p or not q:
        return []
    r = [RF(0) for _ in range(len(p) + len(q) - 1)]
    for i, a in enumerate(p):
        for j, b in enumerate(q):
            r[i + j] = r[i + j] + a * b
    return pnorm(r)


def plin(root):
    root = root if isinstance(root, RF) else RF(Fr(root))
    return [-root, RF(1)]


def pprod(ps):
    r = [RF(1)]
    for p in ps:
        r = pmul(r, p)
    return r


def pev(p, x):
    r = RF(0)
    for c in reversed(p):
        r = r * x + c
    return r


def pdivmod(p, q):
    p = pnorm(p)
    q = pnorm(q)
    p = p[:]
    out = [RF(0) for _ in range(max(len(p) - len(q) + 1, 1))]
    while len(p) >= len(q) and p:
        c = p[-1] / q[-1]
        k = len(p) - len(q)
        out[k] = c
        for i, b in enumerate(q):
            p[i + k] = p[i + k] - c * b
        p = pnorm(p)
    return pnorm(out), p


def pinterp(xs, ys):
    res = []
    for i in range(len(xs)):
        num = [RF(1)]
        den = RF(1)
        for j in range(len(xs)):
            if j != i:
                num = pmul(num, plin(xs[j]))
                den = den * (xs[i] - xs[j])
        res = padd(res, pscal(ys[i] / den, num))
    return res


def det(M):
    M = [row[:] for row in M]
    n = len(M)
    d = RF(1)
    for c in range(n):
        p = None
        for i in range(c, n):
            if not M[i][c].iszero():
                p = i
                break
        if p is None:
            return RF(0)
        if p != c:
            M[c], M[p] = M[p], M[c]
            d = -d
        d = d * M[c][c]
        for i in range(c + 1, n):
            if not M[i][c].iszero():
                f = M[i][c] / M[c][c]
                M[i] = [a - f * b for a, b in zip(M[i], M[c])]
    return d


def sylvester_res(p, q):
    p = pnorm(p)
    q = pnorm(q)
    m = len(p) - 1
    n = len(q) - 1
    N = m + n
    rows = []
    for i in range(n):
        row = [RF(0)] * N
        for k, c in enumerate(reversed(p)):
            row[i + k] = c
        rows.append(row)
    for i in range(m):
        row = [RF(0)] * N
        for k, c in enumerate(reversed(q)):
            row[i + k] = c
        rows.append(row)
    return det(rows)
