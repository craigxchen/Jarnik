#!/usr/bin/env python3
"""Independent exact checker for families.md.

Every normalized constant is computed from exact Gaussian integers:
  * the cluster is divided by its Gaussian gcd (Euclid in Z[i]);
  * equal norms and distinctness are asserted;
  * each point is rotated by a unit towards the first point, and the angle
    atan(Im/Re) of the exact ratio z * conj(z0) is evaluated with Decimal
    (60 digits) by its Taylor series (|Im/Re| < 0.5 in every case below);
  * C = (max angle - min angle) * N^(1/4), N = common squared modulus.
No floating point is used in these values.
"""
import itertools
from decimal import Decimal, getcontext
getcontext().prec = 60

# ---------- Gaussian integer arithmetic ----------
def gmul(a, b): return (a[0]*b[0]-a[1]*b[1], a[0]*b[1]+a[1]*b[0])
def conj(z): return (z[0], -z[1])
def gpow(z, e):
    p = (1, 0)
    for _ in range(e): p = gmul(p, z)
    return p
def ggcd(a, b):
    while b != (0, 0):
        n = b[0]*b[0]+b[1]*b[1]
        t = gmul(a, conj(b))
        q = ((2*t[0]+n)//(2*n), (2*t[1]+n)//(2*n))
        r = gmul(b, q)
        a, b = b, (a[0]-r[0], a[1]-r[1])
    return a
def gdiv(a, b):
    n = b[0]*b[0]+b[1]*b[1]; t = gmul(a, conj(b))
    assert t[0] % n == 0 and t[1] % n == 0
    return (t[0]//n, t[1]//n)

def datan(x):
    # atan for |x| < 0.5 by Taylor series in Decimal
    s = Decimal(0); term = x; k = 0; x2 = x*x
    while True:
        add = term/(2*k+1)
        if k > 0 and abs(add) <= abs(s)*Decimal(10)**-58: break
        if add == 0: break
        s += add if k % 2 == 0 else -add
        term *= x2; k += 1
    return s

def normalized_C(points):
    """Return (C, N, k, gcd) for a cluster of equal-modulus Gaussian integers."""
    g = points[0]
    for p in points[1:]: g = ggcd(g, p)
    P = [gdiv(p, g) for p in points]
    N = P[0][0]**2 + P[0][1]**2
    assert all(p[0]**2+p[1]**2 == N for p in P), "unequal norms"
    assert len(set(P)) == len(P), "points not distinct"
    c0 = conj(P[0]); angs = []
    for p in P:
        q = gmul(p, c0)
        q = max([q, (-q[1], q[0]), (-q[0], -q[1]), (q[1], -q[0])], key=lambda t: t[0])
        x = Decimal(q[1])/Decimal(q[0])
        assert abs(x) < Decimal("0.5")
        angs.append(datan(x))
    # distinct unit classes
    span = max(angs) - min(angs)
    N4 = Decimal(N).sqrt().sqrt()
    return span*N4, N, len(P), g

def best_sub(points, k):
    """least normalized span over k-subsets (consecutive in angle), after the gcd of the
    whole cluster; returns C over the best k consecutive points of the full cluster."""
    g = points[0]
    for p in points[1:]: g = ggcd(g, p)
    P = [gdiv(p, g) for p in points]
    N = P[0][0]**2 + P[0][1]**2
    c0 = conj(P[0]); angs = []
    for p in P:
        q = gmul(p, c0)
        q = max([q, (-q[1], q[0]), (-q[0], -q[1]), (q[1], -q[0])], key=lambda t: t[0])
        angs.append(datan(Decimal(q[1])/Decimal(q[0])))
    angs.sort()
    N4 = Decimal(N).sqrt().sqrt()
    return min(angs[i+k-1]-angs[i] for i in range(len(angs)-k+1))*N4

S5 = Decimal(5).sqrt(); PHI = (1+S5)/2; R2 = Decimal(2).sqrt(); F4 = S5.sqrt()

# ---------- 1. golden (Fibonacci) fixed-shift templates ----------
_f = [0, 1]
def fib(r):
    while len(_f) <= r: _f.append(_f[-1]+_f[-2])
    return _f[r]
def J(r): return (fib(r+1), fib(r))   # norm f_{2r+1}
F = (1, 2)
def golden_rows(n, shifts, widths, parity):
    E = sum(widths)
    rows = [a for a in itertools.product(*[range(w+1) for w in widths]) if sum(a) % 2 == parity]
    es = [(E-2*sum(a))//2 for a in rows]; tot = max(abs(e) for e in es)
    pts = []
    for a, e in zip(rows, es):
        z = gmul(gmul(gpow(F, max(e, 0)), gpow(conj(F), max(-e, 0))), (5**((tot-abs(e))//2), 0))
        for s, w, x in zip(shifts, widths, a):
            b = J(n+s); z = gmul(z, gmul(gpow(b, x), gpow(conj(b), w-x)))
        pts.append(z)
    return pts

def check_golden():
    lim = {8: R2*(4+S5)/F4, 7: R2*PHI**2*F4, 6: R2*PHI**3/F4, 5: R2*F4*PHI}
    print("1. golden templates (blocks J_r = f_{r+1} + i f_r):")
    print("   closed forms: C8=sqrt2(4+sqrt5)/5^(1/4)=%.12f  C7=sqrt2 phi^2 5^(1/4)=%.12f" % (lim[8], lim[7]))
    print("                 C6=sqrt2 phi^3/5^(1/4)=%.12f  C5=sqrt2 5^(1/4) phi=%.12f" % (lim[6], lim[5]))
    for n in (100, 160, 220):   # n = 1 mod 3: two blocks have even norm, gcd norm 4
        assert n % 3 == 1
        pts = golden_rows(n, [0, 1, 2, 3], [1, 1, 1, 1], 1)
        C8, N, k, g = normalized_C(pts)
        C7 = best_sub(pts, 7)
        assert k == 8 and g[0]**2+g[1]**2 == 4
        pts6 = golden_rows(n-1, [0, 1, 2], [1, 2, 1], 1)   # n-1 = 0 mod 3: J_{n} (width 2) has even norm
        C6, N6, k6, g6 = normalized_C(pts6)
        C5 = best_sub(pts6, 5)
        assert k6 == 6
        print("   n=%d: k=8 C=%.12f (err %.1e), k=7 C=%.12f (err %.1e), k=6 C=%.12f (err %.1e), k=5 C=%.12f (err %.1e), log10 R=%.0f"
              % (n, C8, C8-lim[8], C7, C7-lim[7], C6, C6-lim[6], C5, C5-lim[5], N.bit_length()*0.30103/2))
        for kk, CC in ((8, C8), (7, C7), (6, C6), (5, C5)):
            assert abs(CC - lim[kk]) < Decimal("1e-20")
    # small member found by the exhaustive scan: n = 4 gives R^2 = 537606725
    pts = golden_rows(4, [0, 1, 2, 3], [1, 1, 1, 1], 1)
    C8, N, k, g = normalized_C(pts)
    assert N == 537606725 and k == 8
    print("   n=4 member: R^2=%d, 8 points, C=%.6f (the exhaustive scan's k=8 record)" % (N, C8))

# ---------- 2. Cilleruelo-Granville six points after gcd ----------
def check_cg():
    print("2. CG six-point polynomial family prod_j (u-j-i s_j), sum s=0, after Gaussian gcd:")
    for u in (10**6, 10**12+1):
        pts = []
        for S in itertools.combinations(range(4), 2):
            z = (1, 0)
            for j in range(4):
                z = gmul(z, (u-(j+1), -1 if j in S else 1))
            pts.append(z)
        C, N, k, g = normalized_C(pts)
        print("   u=%d: k=%d gcd norm=%d C=%.10f (4 sqrt2 = %.10f)" % (u, k, g[0]**2+g[1]**2, C, 4*R2))
        assert k == 6 and abs(C-4*R2) < Decimal("1e-4")

# ---------- 3. translated Pell eight-point family of the notes ----------
def check_pell8():
    def pell(r):
        x, y = 1, 0
        for _ in range(r): x, y = 9*x+20*y, 4*x+9*y
        return x, y
    def H(r):
        x, y = pell(r); return (x+y, 2*y)
    print("3. translated Pell eight-point family (primitive_eight_point_translated_pell_family.md):")
    for n in (61, 121):
        bl = [H(n-2), H(n), H(n+2), H(n+4)]; pts = []
        for w in itertools.product([0, 1], repeat=4):
            p = sum(w)
            if p % 2 == 0: continue
            z = F if p == 1 else conj(F)
            for b, s in zip(bl, w): z = gmul(z, b if s else conj(b))
            pts.append(z)
        C, N, k, g = normalized_C(pts)
        print("   n=%d: k=%d C=%.6f" % (n, k, C))
        assert k == 8 and Decimal("4.46e7") < C < Decimal("4.47e7")

# ---------- 4. axis six-point family with fixed y0 (Pell conic), C -> 4 ----------
def check_axis6():
    print("4. axis six-point family: n=a^2+y0^2=(a-1)^2+y1^2=(a-2)^2+y2^2, y0=12 fixed:")
    # y2^2 + y0^2 + 2 = 2 y1^2 ; with y0=12: y2^2 - 2 y1^2 = -146 ; a=(y1^2-y0^2+1)/2
    y0 = 12; sols = []
    for y1 in range(1, 3000):
        r = 2*y1*y1 - y0*y0 - 2
        if r > 0:
            y2 = int(Decimal(r).sqrt())
            for t in (y2-1, y2, y2+1):
                if t > 0 and t*t == r and (y1*y1-y0*y0+1) % 2 == 0: sols.append((y1, t))
    assert (1171, 1656) in sols
    # orbit under the unit 3+2sqrt2 acting on y2 + y1 sqrt2
    seeds = [(1171, 1656)]
    out = []
    y1, y2 = 1171, 1656
    for it in range(4):
        a = (y1*y1-y0*y0+1)//2
        assert (y1*y1-y0*y0+1) % 2 == 0 and a > 2
        pts = []
        for d, y in ((0, y0), (1, y1), (2, y2)):
            assert (a-d)**2 + y*y == a*a + y0*y0
            pts += [(a-d, y), (a-d, -y)]
        C, N, k, g = normalized_C(pts)
        print("   a=%d: k=%d C=%.12f" % (a, k, C))
        out.append(C)
        # (y2 + y1 sqrt2)(3+2 sqrt2) = (3 y2 + 4 y1) + (2 y2 + 3 y1) sqrt2 ; parity of y1 preserved?
        y2, y1 = 3*y2 + 4*y1, 2*y2 + 3*y1
        # restore odd y1 by applying the unit twice if needed
        if (y1*y1 - y0*y0 + 1) % 2: y2, y1 = 3*y2 + 4*y1, 2*y2 + 3*y1
    assert all(c > 4 for c in out) and out[-1] - 4 < Decimal("1e-8")

# ---------- 5. no circle has lattice points at four consecutive abscissae ----------
def check_mod8():
    sq = {x*x % 8 for x in range(8)}
    bad = [(n, a) for n in range(8) for a in range(8)
           if all((n-(a-j)**2) % 8 in sq for j in range(4))]
    print("5. residues (n,a) mod 8 with n-(a-j)^2 a square mod 8 for j=0..3:", bad)
    assert bad == []

# ---------- 6. sporadic record clusters from the scans ----------
def points_from(n, xs):
    pts = []
    for x in xs:
        r = n - x*x; y = int(Decimal(r).sqrt())
        y = next(t for t in (y-1, y, y+1) if t >= 0 and t*t == r)
        pts.append((x, y))
    return pts
def check_sporadic():
    print("6. record clusters (exact):")
    cases = [
        (32988780417125, [(5743582, 6799), (5743583, 5894), (5743585, 3430), (5743586, 527),
                          (5743586, -527), (5743585, -3430), (5743583, -5894), (5743582, -6799)], "k=8 axis D=(0,1,3,4)"),
        (348984761607530165, [(590749322-d, sg*y) for d, y in ((0, 12809), (1, 36682), (3, 60898), (4, 69929)) for sg in (1, -1)], "k=8 axis D=(0,1,3,4), second non-golden point"),
        (235370525, [(14650, 4555), (14629, 4622), (14555, 4850), (14554, 4853), (14491, 5038),
                     (14453, 5146), (14443, 5174)], "k=7 sporadic"),
        (469977431545, [(685547, 1656), (685548, 1171), (685549, 12), (685549, -12), (685548, -1171),
                        (685547, -1656)], "k=6 axis D=(0,1,2)"),
        (683617125425, [(706952, 428761), (706544, 429433), (706513, 429484), (706087, 430184),
                        (706009, 430312)], "k=5 sporadic"),
        (276108509, None, "k=3"),
    ]
    for n, pts, label in cases:
        if pts is None: continue
        assert all(x*x+y*y == n for x, y in pts)
        C, N, k, g = normalized_C(pts)
        print("   R^2=%d %s: primitive gcd norm %d, k=%d, C=%.9f" % (n, label, g[0]**2+g[1]**2, k, C))

if __name__ == "__main__":
    check_golden(); check_cg(); check_pell8(); check_axis6(); check_mod8(); check_sporadic()
    print("ALL CHECKS PASSED")
