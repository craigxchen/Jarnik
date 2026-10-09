"""Exact univariate and multivariate polynomial helpers over Q (fractions.Fraction).

Univariate polynomials: lists of Fractions, low degree first, trailing zeros stripped.
Multivariate polynomials: dict {exponent tuple: Fraction}.
No floating point is used in any exact routine; `rational_roots` uses floats only to
propose candidates, every candidate is verified exactly.
"""
from fractions import Fraction as Fr
from math import gcd
import itertools

# ----------------------------------------------------------------- univariate

def norm(p):
    p = [Fr(c) for c in p]
    while p and p[-1] == 0:
        p.pop()
    return p

def deg(p):
    return len(norm(p)) - 1

def add(p, q):
    n = max(len(p), len(q))
    return norm([(p[i] if i < len(p) else 0) + (q[i] if i < len(q) else 0) for i in range(n)])

def sub(p, q):
    return add(p, [-c for c in q])

def mul(p, q):
    if not p or not q:
        return []
    r = [Fr(0)] * (len(p) + len(q) - 1)
    for i, a in enumerate(p):
        if a == 0:
            continue
        for j, b in enumerate(q):
            r[i + j] += a * b
    return norm(r)

def scal(c, p):
    return norm([c * a for a in p])

def ev(p, x):
    r = Fr(0)
    for c in reversed(p):
        r = r * x + c
    return r

def divmod_(p, q):
    p = norm(p); q = norm(q)
    if not q:
        raise ZeroDivisionError
    quo = [Fr(0)] * max(len(p) - len(q) + 1, 1)
    r = p[:]
    while len(r) >= len(q) and r:
        c = r[-1] / q[-1]
        k = len(r) - len(q)
        quo[k] = c
        for i, b in enumerate(q):
            r[i + k] -= c * b
        r = norm(r)
    return norm(quo), r

def pgcd(p, q):
    p = norm(p); q = norm(q)
    while q:
        _, r = divmod_(p, q)
        p, q = q, r
    if p:
        p = scal(1 / p[-1], p)
    return p

def deriv(p):
    return norm([i * p[i] for i in range(1, len(p))])

def linpoly(a, b):
    """a*x + b"""
    return norm([Fr(b), Fr(a)])

def prod(ps):
    r = [Fr(1)]
    for p in ps:
        r = mul(r, p)
    return r

def primitive_int(p):
    """Scale p to a primitive integer polynomial (list of ints)."""
    p = norm(p)
    if not p:
        return []
    den = 1
    for c in p:
        den = den * c.denominator // gcd(den, c.denominator)
    ints = [int(c * den) for c in p]
    g = 0
    for c in ints:
        g = gcd(g, c)
    return [c // g for c in ints]

def rational_roots(p, maxden=10**12):
    """All rational roots of p (exact). Candidates from numpy roots + continued fractions."""
    import numpy as np
    p = norm(p)
    if len(p) <= 1:
        return []
    roots = set()
    # strip zero roots
    while p and p[0] == 0:
        roots.add(Fr(0))
        p = p[1:]
    if len(p) <= 1:
        return sorted(roots)
    # square-free part for numerical stability
    g = pgcd(p, deriv(p))
    sf, _ = divmod_(p, g) if len(g) > 1 else (p, [])
    co = [float(c) for c in reversed(sf)]
    try:
        cand = np.roots(co)
    except Exception:
        cand = []
    for z in cand:
        if abs(z.imag) > 1e-6 * max(1.0, abs(z.real)):
            continue
        x = Fr(z.real).limit_denominator(maxden)
        for y in (x, Fr(z.real).limit_denominator(10**6), Fr(z.real).limit_denominator(10**3)):
            if ev(sf, y) == 0:
                roots.add(y)
    return sorted(roots)

def poly_str(p, var='x'):
    p = norm(p)
    if not p:
        return '0'
    terms = []
    for i, c in enumerate(p):
        if c == 0:
            continue
        terms.append(f"({c})*{var}^{i}" if i else f"({c})")
    return ' + '.join(terms)

# --------------------------------------------------------------- multivariate

class MP:
    """Multivariate polynomial over Q, nvars fixed."""
    __slots__ = ('n', 't')

    def __init__(self, n, t=None):
        self.n = n
        self.t = {} if t is None else {k: Fr(v) for k, v in t.items() if v != 0}

    @staticmethod
    def const(n, c):
        return MP(n, {(0,) * n: Fr(c)} if c != 0 else {})

    @staticmethod
    def var(n, i):
        e = [0] * n; e[i] = 1
        return MP(n, {tuple(e): Fr(1)})

    def __add__(self, o):
        if not isinstance(o, MP):
            o = MP.const(self.n, o)
        t = dict(self.t)
        for k, v in o.t.items():
            t[k] = t.get(k, 0) + v
            if t[k] == 0:
                del t[k]
        return MP(self.n, t)
    __radd__ = __add__

    def __neg__(self):
        return MP(self.n, {k: -v for k, v in self.t.items()})

    def __sub__(self, o):
        return self + (-o if isinstance(o, MP) else MP.const(self.n, -o))

    def __rsub__(self, o):
        return (-self) + o

    def __mul__(self, o):
        if not isinstance(o, MP):
            o = Fr(o)
            if o == 0:
                return MP(self.n)
            return MP(self.n, {k: v * o for k, v in self.t.items()})
        t = {}
        for k1, v1 in self.t.items():
            for k2, v2 in o.t.items():
                k = tuple(a + b for a, b in zip(k1, k2))
                t[k] = t.get(k, 0) + v1 * v2
        return MP(self.n, {k: v for k, v in t.items() if v != 0})
    __rmul__ = __mul__

    def __pow__(self, e):
        r = MP.const(self.n, 1)
        for _ in range(e):
            r = r * self
        return r

    def iszero(self):
        return not self.t

    def totdeg(self):
        return max((sum(k) for k in self.t), default=-1)

    def deg_in(self, i):
        return max((k[i] for k in self.t), default=-1)

    def subs(self, vals):
        """vals: dict index -> Fraction; returns MP in same n (substituted vars have exponent 0)."""
        t = {}
        for k, v in self.t.items():
            c = v
            kk = list(k)
            for i, x in vals.items():
                if kk[i]:
                    c *= Fr(x) ** kk[i]
                    kk[i] = 0
            kk = tuple(kk)
            t[kk] = t.get(kk, 0) + c
        return MP(self.n, {k: v for k, v in t.items() if v != 0})

    def evalf(self, vals):
        r = Fr(0)
        for k, v in self.t.items():
            c = v
            for i, e in enumerate(k):
                if e:
                    c *= Fr(vals[i]) ** e
            r += c
        return r

    def as_univariate(self, i):
        """Assuming only variable i remains, return univariate list."""
        d = self.deg_in(i)
        out = [Fr(0)] * (d + 1 if d >= 0 else 0)
        for k, v in self.t.items():
            for j, e in enumerate(k):
                if j != i and e:
                    raise ValueError('other variables present')
            out[k[i]] += v
        return norm(out)


def det_poly(M, mul_=mul, add_=add, sub_=sub):
    """Determinant of a square matrix of univariate polynomials via fraction-free expansion by
    permutations (fine up to 6x6)."""
    n = len(M)
    total = []
    for perm in itertools.permutations(range(n)):
        # sign
        sgn = 1
        p = list(perm)
        for i in range(n):
            while p[i] != i:
                j = p[i]
                p[i], p[j] = p[j], p[i]
                sgn = -sgn
        term = [Fr(1)]
        zero = False
        for i in range(n):
            e = M[i][perm[i]]
            if not norm(e):
                zero = True
                break
            term = mul_(term, e)
        if zero:
            continue
        total = add_(total, term) if sgn > 0 else sub_(total, term)
    return total


def det_frac(M):
    """Exact determinant of a square matrix of Fractions (Gaussian elimination)."""
    A = [[Fr(x) for x in row] for row in M]
    n = len(A)
    d = Fr(1)
    for c in range(n):
        piv = None
        for r in range(c, n):
            if A[r][c] != 0:
                piv = r
                break
        if piv is None:
            return Fr(0)
        if piv != c:
            A[c], A[piv] = A[piv], A[c]
            d = -d
        d *= A[c][c]
        for r in range(c + 1, n):
            f = A[r][c] / A[c][c]
            if f:
                for k in range(c, n):
                    A[r][k] -= f * A[c][k]
    return d


def nullspace_frac(M):
    """Exact right nullspace basis of a Fraction matrix (list of rows)."""
    A = [[Fr(x) for x in row] for row in M]
    m = len(A); n = len(A[0]) if m else 0
    pivcols = []
    r = 0
    for c in range(n):
        piv = None
        for i in range(r, m):
            if A[i][c] != 0:
                piv = i
                break
        if piv is None:
            continue
        A[r], A[piv] = A[piv], A[r]
        pv = A[r][c]
        A[r] = [x / pv for x in A[r]]
        for i in range(m):
            if i != r and A[i][c] != 0:
                f = A[i][c]
                A[i] = [a - f * b for a, b in zip(A[i], A[r])]
        pivcols.append(c)
        r += 1
        if r == m:
            break
    free = [c for c in range(n) if c not in pivcols]
    basis = []
    for f in free:
        v = [Fr(0)] * n
        v[f] = Fr(1)
        for i, pc in enumerate(pivcols):
            v[pc] = -A[i][f]
        basis.append(v)
    return basis


def rank_frac(M):
    if not M:
        return 0
    return len(M[0]) - len(nullspace_frac(M))


def rank_mod(M, p):
    """Rank of an integer/Fraction matrix modulo prime p (entries must be p-integral)."""
    A = []
    for row in M:
        rr = []
        for x in row:
            x = Fr(x)
            if x.denominator % p == 0:
                raise ValueError('not p-integral')
            rr.append(x.numerator * pow(x.denominator, -1, p) % p)
        A.append(rr)
    m = len(A); n = len(A[0]) if m else 0
    r = 0
    for c in range(n):
        piv = None
        for i in range(r, m):
            if A[i][c]:
                piv = i
                break
        if piv is None:
            continue
        A[r], A[piv] = A[piv], A[r]
        inv = pow(A[r][c], -1, p)
        A[r] = [x * inv % p for x in A[r]]
        for i in range(m):
            if i != r and A[i][c]:
                f = A[i][c]
                A[i] = [(a - f * b) % p for a, b in zip(A[i], A[r])]
        r += 1
        if r == m:
            break
    return r
