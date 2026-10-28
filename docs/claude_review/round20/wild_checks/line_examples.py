# Heights of a point u in P^n(Q(i)) and of its Q-line ell_u = span(Re u, Im u)  (wild.md Thm B, Prop C).
# h(u)   = log max|u_x|  after dividing by the Gaussian gcd of the coordinates (exact gcd in Z[i]).
# h(ell) = log max |x_a y_b - x_b y_a| - log gcd(all minors)    (Weil height of the Pluecker point over Q).
# All arithmetic is exact (Python ints); logs are taken at the end.
import math, sys
from itertools import combinations
from math import gcd

def gmod(a, b):  # Gaussian remainder, a,b tuples
    (x, y), (u, v) = a, b
    n = u*u + v*v
    # a * conj(b) / n rounded
    re = x*u + y*v; im = y*u - x*v
    qr = (2*re + n) // (2*n); qi = (2*im + n) // (2*n)
    return (x - (qr*u - qi*v), y - (qr*v + qi*u))

def ggcd(a, b):
    while b != (0, 0):
        a, b = b, gmod(a, b)
    return a

def gdiv(a, b):
    (x, y), (u, v) = a, b
    n = u*u + v*v
    re = x*u + y*v; im = y*u - x*v
    assert re % n == 0 and im % n == 0
    return (re//n, im//n)

def gmul(a, b):
    return (a[0]*b[0] - a[1]*b[1], a[0]*b[1] + a[1]*b[0])

def gpow(a, k):
    r = (1, 0)
    for _ in range(k): r = gmul(r, a)
    return r

def heights(pts):
    g = (0, 0)
    for p in pts: g = ggcd(g, p)
    pts = [gdiv(p, g) for p in pts]
    hu = max(0.5*math.log(x*x + y*y) for x, y in pts)
    minors = [pts[a][0]*pts[b][1] - pts[b][0]*pts[a][1] for a, b in combinations(range(len(pts)), 2)]
    G = 0
    for m in minors: G = gcd(G, m)
    mx = max(abs(m) for m in minors)
    hl = math.log(mx) - math.log(G)
    return hu, hl, G

if __name__ == "__main__":
    print("== Prop C: rank-one points u=(1,g,g^b), g=A+i (subgroup of codim n-1 in G_m^n) ==")
    for b in (2, 3, 5):
        for A in (10, 10**3, 10**6, 10**12):
            g = (A, 1)
            pts = [gpow(g, k) for k in ([0, 1, b] if b > 1 else [0, 1])]
            hu, hl, G = heights(pts)
            print(f"b={b} A={A:>14}: h(u)={hu:10.4f} h(ell)={hl:10.4f} ratio={hu/hl:.6f} gcd(minors)={G}")
    for A in (10**3, 10**6):
        g = (A, 1); n = 8
        pts = [gpow(g, k) for k in range(n+1)]
        hu, hl, G = heights(pts)
        print(f"u=(g^0..g^{n}), A={A}: h(u)={hu:.4f} h(ell)={hl:.4f} ratio={hu/hl:.6f} G={G}")
