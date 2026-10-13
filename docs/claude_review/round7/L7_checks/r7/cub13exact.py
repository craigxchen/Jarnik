"""EXACT construction of cub13 curves in the Kapranov psi_7 model (works over any field via
number-like objects: Fraction, nf.KE, nf.AD).

Model: P^4 with p1..p5 = e0..e4, p6 = (1,...,1).  A point (c1:...:c5) is the configuration
(y1,...,y7) = (c1,...,c5, 0, inf) of 7 points on P^1 (labels 6 -> 0, 7 -> inf).
Curve: c(t) = (h, X2, X3, X4, X5), cubics.  Gauge: t_126 = inf (lam3=lam4=lam5=1), r45 = 0, r24 = 1.
Free parameters p = (r25, r34, r35, w); w = t_456.  Membership in the cub13 family: residual(p) = 0."""

def pmul(a, b):
    r = [0] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b): r[i + j] = r[i + j] + x * y
    return r
def padd(a, b):
    n = max(len(a), len(b))
    return [(a[i] if i < len(a) else 0) + (b[i] if i < len(b) else 0) for i in range(n)]
def pscale(c, a): return [c * x for x in a]
def lin(r): return [-r, 1]
def ev(a, x):
    r = 0
    for c in reversed(a): r = r * x + c
    return r
def prod_lin(rs, lead=1):
    p = [lead]
    for r in rs: p = pmul(p, lin(r))
    return p

def build(r25, r34, r35, w, one=1):
    """returns dict with polynomials X2..X5, h (low-first coefficient lists), contact parameters,
    and residual (zero iff the curve meets plane 456 at t = w)."""
    r24 = one * 1; r45 = one * 0
    v = (r25 * r35 - r24 * r34) / (r25 + r35 - r24 - r34)
    r23 = ((r25 * r35 - r34) * (r25 * r34 * r35 - r25 * r34 - 2 * r25 * r35 * w + r25 * r35 + r25 * w * w
            - r34 * r35 - r34 * w * w + 2 * r34 * w + r35 * w * w - w * w)
           / ((r25 - 1) * (r34 - w) * (r35 - w) * (r25 - r34 + r35 - 1)))
    X4 = prod_lin([r24, r34, r45], one)
    X5 = prod_lin([r25, r35, r45], one)
    X3 = prod_lin([r23, r34, r35], one)
    lam2 = ev(X4, v) / ((v - r23) * (v - r24) * (v - r25))
    X2 = prod_lin([r23, r24, r25], lam2)
    def other_root(P, known, deg):
        # P has degree deg (coefficient list possibly with a vanishing top coefficient beyond deg)
        s = -P[deg - 1] / P[deg]
        for k in known: s = s - k
        return s
    D35 = padd(X3, pscale(-1, X5)); D34 = padd(X3, pscale(-1, X4))
    D25 = padd(X2, pscale(-1, X5)); D24 = padd(X2, pscale(-1, X4))
    t246 = other_root(D35, [r35], 2)
    t256 = other_root(D34, [r34], 2)
    t346 = other_root(D25, [r25, v], 3)
    t356 = other_root(D24, [r24, v], 3)
    pts = [(t246, ev(X3, t246)), (t256, ev(X3, t256)), (t346, ev(X2, t346)), (t356, ev(X2, t356))]
    # Lagrange interpolation of the cubic h through the 4 points
    h = [0]
    for i, (xi, yi) in enumerate(pts):
        num = [one * 1]; den = one * 1
        for j, (xj, _) in enumerate(pts):
            if j != i:
                num = pmul(num, lin(xj)); den = den * (xi - xj)
        h = padd(h, pscale(yi / den, num))
    residual = ev(h, w) - ev(X2, w)
    return dict(X2=X2, X3=X3, X4=X4, X5=X5, h=h, residual=residual, v=v, r23=r23, lam2=lam2,
                contacts=dict(t246=t246, t256=t256, t346=t346, t356=t356, t456=w, t126='inf',
                              t123=r45, t135=r24), r=dict(r23=r23, r24=r24, r25=r25, r34=r34, r35=r35, r45=r45))

def coords(d):
    """the five cubics c1..c5 = (h, X2, X3, X4, X5) <-> labels 1..5"""
    return [d['h'], d['X2'], d['X3'], d['X4'], d['X5']]
