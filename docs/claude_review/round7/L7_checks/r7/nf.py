"""Exact arithmetic in a number field K = Q[x]/(m(x)) (m irreducible over Q), plus forward-mode
automatic differentiation over K and univariate polynomials over K (gcd, derivative).
Pure Python + fractions; no floating point anywhere."""
from fractions import Fraction as Fr

def _trim(a):
    a = list(a)
    while a and a[-1] == 0:
        a.pop()
    return a

def _pmul(a, b):
    if not a or not b: return []
    r = [Fr(0)] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        if x == 0: continue
        for j, y in enumerate(b):
            r[i + j] += x * y
    return _trim(r)

def _pdivmod(a, b):
    a = _trim(a); b = _trim(b)
    q = [Fr(0)] * max(len(a) - len(b) + 1, 1)
    a = a[:]
    while len(a) >= len(b) and a:
        c = a[-1] / b[-1]; k = len(a) - len(b); q[k] = c
        for i, y in enumerate(b): a[i + k] -= c * y
        a = _trim(a)
    return _trim(q), a

class Field:
    def __init__(self, m_low_first):
        m = [Fr(c) for c in m_low_first]
        self.m = [c / m[-1] for c in m]  # monic
        self.n = len(self.m) - 1
    def __call__(self, v):
        return KE(self, v)
    def gen(self):
        return KE(self, [Fr(0), Fr(1)])

class KE:
    __slots__ = ('F', 'c')
    def __init__(self, F, v):
        self.F = F
        if isinstance(v, KE): self.c = v.c[:]
        elif isinstance(v, (list, tuple)): self.c = _pdivmod([Fr(x) for x in v], F.m)[1] if len(v) > F.n else _trim([Fr(x) for x in v])
        else: self.c = _trim([Fr(v)])
    def _co(self, o):
        return o if isinstance(o, KE) else KE(self.F, o)
    def __add__(s, o):
        o = s._co(o); n = max(len(s.c), len(o.c))
        return KE(s.F, [(s.c[i] if i < len(s.c) else 0) + (o.c[i] if i < len(o.c) else 0) for i in range(n)])
    __radd__ = __add__
    def __neg__(s): return KE(s.F, [-x for x in s.c])
    def __sub__(s, o): return s + (-s._co(o))
    def __rsub__(s, o): return s._co(o) - s
    def __mul__(s, o):
        o = s._co(o)
        return KE(s.F, _pdivmod(_pmul(s.c, o.c), s.F.m)[1])
    __rmul__ = __mul__
    def inv(s):
        if not s.c: raise ZeroDivisionError('zero in K')
        # extended Euclid over Q[x]: track u with u*s == r (mod m)
        r0, r1 = s.F.m[:], s.c[:]
        u0, u1 = [], [Fr(1)]
        while len(r1) > 1:
            q, r = _pdivmod(r0, r1)
            r0, r1 = r1, r
            u0, u1 = u1, _psub(u0, _pmul(q, u1))
            if not r1: raise ZeroDivisionError('not invertible (m reducible?)')
        c = r1[0]
        return KE(s.F, [x / c for x in u1])
    def __truediv__(s, o): return s * s._co(o).inv()
    def __rtruediv__(s, o): return s._co(o) * s.inv()
    def __pow__(s, k):
        r = KE(s.F, 1)
        for _ in range(k): r = r * s
        return r
    def __eq__(s, o): return (s - s._co(o)).c == []
    def iszero(s): return s.c == []
    def __repr__(s): return 'K(%s)' % [str(x) for x in s.c]

def _psub(a, b):
    n = max(len(a), len(b))
    return _trim([(a[i] if i < len(a) else 0) - (b[i] if i < len(b) else 0) for i in range(n)])

class AD:
    """value in K, gradient = list of K (length nv)."""
    __slots__ = ('v', 'g')
    def __init__(self, v, g): self.v, self.g = v, g
    @staticmethod
    def const(F, x, nv): return AD(KE(F, x) if not isinstance(x, KE) else x, [KE(F, 0)] * nv)
    def _co(s, o): return o if isinstance(o, AD) else AD(s.v._co(o), [KE(s.v.F, 0)] * len(s.g))
    def __add__(s, o): o = s._co(o); return AD(s.v + o.v, [a + b for a, b in zip(s.g, o.g)])
    __radd__ = __add__
    def __neg__(s): return AD(-s.v, [-a for a in s.g])
    def __sub__(s, o): return s + (-s._co(o))
    def __rsub__(s, o): return s._co(o) - s
    def __mul__(s, o): o = s._co(o); return AD(s.v * o.v, [s.v * b + o.v * a for a, b in zip(s.g, o.g)])
    __rmul__ = __mul__
    def __truediv__(s, o):
        o = s._co(o); iv = o.v.inv(); q = s.v * iv
        return AD(q, [(a - q * b) * iv for a, b in zip(s.g, o.g)])
    def __rtruediv__(s, o): return s._co(o) / s

# ---- univariate polynomials over K: lists of KE, low first --------------------------------
def ptrim(a):
    a = list(a)
    while a and a[-1].iszero(): a.pop()
    return a
def padd(a, b):
    n = max(len(a), len(b)); F = (a or b)[0].F
    return ptrim([(a[i] if i < len(a) else KE(F, 0)) + (b[i] if i < len(b) else KE(F, 0)) for i in range(n)])
def pneg(a): return [-x for x in a]
def psub(a, b): return padd(a, pneg(b))
def pmul(a, b):
    if not a or not b: return []
    F = a[0].F; r = [KE(F, 0)] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b): r[i + j] = r[i + j] + x * y
    return ptrim(r)
def pdivmod(a, b):
    a = ptrim(a); b = ptrim(b); F = b[0].F
    q = [KE(F, 0)] * max(len(a) - len(b) + 1, 1); a = a[:]
    lead_inv = b[-1].inv()
    while len(a) >= len(b) and a:
        c = a[-1] * lead_inv; k = len(a) - len(b); q[k] = c
        for i, y in enumerate(b): a[i + k] = a[i + k] - c * y
        a = ptrim(a)
    return ptrim(q), a
def pgcd(a, b):
    a, b = ptrim(a), ptrim(b)
    while b:
        a, b = b, pdivmod(a, b)[1]
    if not a: return a
    li = a[-1].inv()
    return [x * li for x in a]
def pderiv(a):
    return ptrim([a[i] * i for i in range(1, len(a))])
def pdeg(a): return len(ptrim(a)) - 1
def peval(a, x):
    r = KE(x.F, 0) if isinstance(x, KE) else 0
    for c in reversed(a): r = r * x + c
    return r
