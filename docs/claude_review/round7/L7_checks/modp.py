"""Arithmetic in F_p (p = 2^61 - 1) for fast exploration of degrees (not used as proof)."""
PR = (1 << 61) - 1
class M:
    __slots__ = ('v',)
    def __init__(s, v):
        if isinstance(v, M): s.v = v.v
        else:
            try:
                from fractions import Fraction
                if isinstance(v, Fraction):
                    s.v = v.numerator % PR * pow(v.denominator % PR, PR - 2, PR) % PR; return
            except Exception: pass
            s.v = int(v) % PR
    def __add__(s, o): return M(s.v + M(o).v)
    __radd__ = __add__
    def __sub__(s, o): return M(s.v - M(o).v)
    def __rsub__(s, o): return M(M(o).v - s.v)
    def __mul__(s, o): return M(s.v * M(o).v)
    __rmul__ = __mul__
    def __neg__(s): return M(-s.v)
    def __truediv__(s, o):
        o = M(o)
        if o.v == 0: raise ZeroDivisionError
        return M(s.v * pow(o.v, PR - 2, PR))
    def __rtruediv__(s, o): return M(o) / s
    def __pow__(s, k): return M(pow(s.v, k, PR))
    def __eq__(s, o): return s.v == M(o).v
    def __hash__(s): return hash(s.v)
    def __repr__(s): return 'M(%d)' % s.v

def nullspace(Mx):
    Mx = [[M(a) for a in row] for row in Mx]
    nr, nc = len(Mx), len(Mx[0]); piv = []; r = 0
    for c in range(nc):
        p = None
        for i in range(r, nr):
            if Mx[i][c].v != 0: p = i; break
        if p is None: continue
        Mx[r], Mx[p] = Mx[p], Mx[r]
        inv = M(1) / Mx[r][c]
        Mx[r] = [x * inv for x in Mx[r]]
        for i in range(nr):
            if i != r and Mx[i][c].v != 0:
                f = Mx[i][c]
                Mx[i] = [a - f * b for a, b in zip(Mx[i], Mx[r])]
        piv.append(c); r += 1
        if r == nr: break
    free = [c for c in range(nc) if c not in piv]
    basis = []
    for fcol in free:
        v = [M(0)] * nc; v[fcol] = M(1)
        for i, c in enumerate(piv): v[c] = -Mx[i][fcol]
        basis.append(v)
    return basis

def pnorm(p):
    p = list(p)
    while p and p[-1].v == 0: p.pop()
    return p
def pdivmod(p, q):
    p = pnorm(p)[:]; q = pnorm(q)
    out = [M(0)] * max(len(p) - len(q) + 1, 1)
    while len(p) >= len(q) and p:
        c = p[-1] / q[-1]; k = len(p) - len(q); out[k] = c
        for i, b in enumerate(q): p[i + k] = p[i + k] - c * b
        p = pnorm(p)
    return pnorm(out), p
def pgcd(p, q):
    p, q = pnorm(p), pnorm(q)
    while q:
        _, r = pdivmod(p, q); p, q = q, r
    inv = M(1) / p[-1]
    return [c * inv for c in p]
def pev(p, x):
    r = M(0)
    for c in reversed(p): r = r * x + c
    return r

def ratdeg(f, D=40, seed=12345):
    """degrees (deg num, deg den) of the rational function f over F_p, plus num, den."""
    import random
    rnd = random.Random(seed)
    pts = []
    while len(pts) < 2 * D + 6:
        x = M(rnd.randrange(PR))
        try: pts.append((x, M(f(x))))
        except ZeroDivisionError: pass
    rows = [[x ** i for i in range(D + 1)] + [-(y * x ** j) for j in range(D + 1)] for x, y in pts]
    ns = nullspace(rows)
    if not ns: return None
    v = ns[0]
    p = pnorm(v[:D + 1]); q = pnorm(v[D + 1:])
    g = pgcd(p, q)
    p, _ = pdivmod(p, g); q, _ = pdivmod(q, g)
    # check
    for _ in range(5):
        x = M(rnd.randrange(PR))
        try:
            assert pev(p, x) == f(x) * pev(q, x)
        except ZeroDivisionError: pass
    return len(p) - 1, len(q) - 1, p, q

def roots_in_set(p, cands):
    return [c for c in cands if pev(p, M(c)).v == 0]
