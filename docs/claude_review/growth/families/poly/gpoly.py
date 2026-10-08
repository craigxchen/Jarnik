# Exact helpers for polynomial templates z(u) = prod (u - beta_m), beta in Z[i] (as (x,y) int pairs).
from fractions import Fraction
import math

def gmul(a, b):
    return (a[0]*b[0]-a[1]*b[1], a[0]*b[1]+a[1]*b[0])

def gpow_im(b, r):
    p = (1, 0)
    for _ in range(r):
        p = gmul(p, b)
    return p[1]

def poly_from_roots(roots):
    # coefficients (list of Gaussian ints, low->high) of prod (u - beta)
    c = [(1, 0)]
    for (x, y) in roots:
        n = [(0, 0)]*(len(c)+1)
        for i, a in enumerate(c):
            # a*u
            n[i+1] = (n[i+1][0]+a[0], n[i+1][1]+a[1])
            # -beta*a
            t = gmul(a, (x, y))
            n[i] = (n[i][0]-t[0], n[i][1]-t[1])
        c = n
    return c

def peval(c, u):
    re = 0; im = 0
    for a in reversed(c):
        re, im = re*u + a[0], im*u + a[1]
    return (re, im)

def ggcd(a, b):
    # Euclid in Z[i]
    while b != (0, 0):
        n = b[0]*b[0]+b[1]*b[1]
        # a/b = a*conj(b)/n
        t = gmul(a, (b[0], -b[1]))
        qx = (2*t[0] + n)//(2*n); qy = (2*t[1] + n)//(2*n)
        r = gmul(b, (qx, qy))
        a, b = b, (a[0]-r[0], a[1]-r[1])
    return a

def gdiv(a, b):
    n = b[0]*b[0]+b[1]*b[1]
    t = gmul(a, (b[0], -b[1]))
    assert t[0] % n == 0 and t[1] % n == 0
    return (t[0]//n, t[1]//n)

def cluster_C(points, return_parts=False):
    """Normalized arc constant C = (angular span) * sqrt(R) of a cluster of Gaussian
    integers of one modulus, after dividing by their Gaussian gcd.  Each point is
    rotated by a unit to lie closest to the first point; valid for clusters whose
    span is < pi/4 (always the case here).  Ratios Im/Re are computed from exact
    integers, so there is no cancellation error."""
    g = points[0]
    for p in points[1:]:
        g = ggcd(g, p)
    pts = [gdiv(p, g) for p in points]
    N = pts[0][0]**2 + pts[0][1]**2
    assert all(p[0]**2+p[1]**2 == N for p in pts)
    p0c = (pts[0][0], -pts[0][1])
    angs = []
    for p in pts:
        q = gmul(p, p0c)
        # rotate by unit so that Re q is maximal
        cands = [q, (-q[1], q[0]), (-q[0], -q[1]), (q[1], -q[0])]
        q = max(cands, key=lambda t: t[0])
        angs.append(math.atan(q[1]/q[0]) if q[1] != 0 else 0.0)
    span = max(angs) - min(angs)
    # sqrt(R) = N^(1/4); use integer sqrt for large N
    r4 = math.isqrt(math.isqrt(N))
    sq = float(r4) if r4 > 10**8 else N**0.25
    C = span*sq
    if return_parts:
        return C, N, len(set(pts)), g
    return C, N, len(set(pts))
