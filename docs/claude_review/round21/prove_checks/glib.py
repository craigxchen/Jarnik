"""Small exact Gaussian-integer helpers (Python ints only)."""
from sympy import factorint, primerange
import mpmath
mpmath.mp.dps = 50

def gmul(a, b):
    return (a[0]*b[0]-a[1]*b[1], a[0]*b[1]+a[1]*b[0])

def gconj(a):
    return (a[0], -a[1])

def gnorm(a):
    return a[0]*a[0]+a[1]*a[1]

def gdivmod_round(a, b):
    # nearest-integer quotient a/b
    n = gnorm(b)
    num = gmul(a, gconj(b))
    def rnd(x):  # round x/n to nearest
        return (2*x + n)//(2*n)
    q = (rnd(num[0]), rnd(num[1]))
    qb = gmul(q, b)
    return q, (a[0]-qb[0], a[1]-qb[1])

def ggcd(a, b):
    while b != (0, 0):
        _, r = gdivmod_round(a, b)
        a, b = b, r
    return a

def gprime(p):
    """a+bi with a>b>0, a^2+b^2=p (p = 1 mod 4)."""
    a = 1
    while True:
        b2 = p - a*a
        b = int(round(b2 ** 0.5))
        while b*b > b2: b -= 1
        while (b+1)*(b+1) <= b2: b += 1
        if b*b == b2 and b > 0:
            return (max(a, b), min(a, b))
        a += 1

def points_mod_units(n):
    """All Gaussian integers of norm n with angle in [0, pi/2), for n with only split prime factors."""
    f = factorint(n)
    pts = [(1, 0)]
    for p, e in f.items():
        assert p % 4 == 1, n
        pi = gprime(p); pib = gconj(pi)
        pw = [(1, 0)]
        for _ in range(e): pw.append(gmul(pw[-1], pi))
        pwb = [(1, 0)]
        for _ in range(e): pwb.append(gmul(pwb[-1], pib))
        new = []
        for z in pts:
            for a in range(e+1):
                new.append(gmul(z, gmul(pw[a], pwb[e-a])))
        pts = new
    out = []
    for z in pts:
        x, y = z
        for _ in range(4):  # rotate into first quadrant: x > 0, y >= 0
            if x > 0 and y >= 0: break
            x, y = -y, x
        out.append((x, y))
    out.sort(key=lambda z: mpmath.atan2(z[1], z[0]))
    return out, f

def windows(n, kmin, Cmax):
    """Yield (k, C, points) for every cyclic window of k consecutive points with chord constant <= Cmax.
    C = max pairwise chord / n^(1/4) (exact squared chords, mpmath root)."""
    pts, f = points_mod_units(n)
    d = len(pts)
    q = mpmath.mpf(n) ** mpmath.mpf('0.25')
    res = []
    for i in range(d):
        for k in range(kmin, d+1):
            W = [pts[(i+t) % d] for t in range(k)]
            # rotate wrapped points by i to keep them adjacent (angles mod pi/2): multiply by i = (0,1) maps
            # first-quadrant point to second quadrant; to compare chords across the wrap use z*i for the
            # points that wrapped: they sit just beyond angle pi/2.
            Wr = [W[t] if i+t < d else (-W[t][1], W[t][0]) for t in range(k)]
            m = 0
            for a in range(k):
                for b in range(a+1, k):
                    dx = Wr[a][0]-Wr[b][0]; dy = Wr[a][1]-Wr[b][1]
                    m = max(m, dx*dx+dy*dy)
            C = mpmath.sqrt(m)/q
            if C > Cmax: break
            res.append((k, C, Wr))
    return res, f
