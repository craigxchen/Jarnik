import math, sys
from gpoly import gmul, ggcd, gdiv
def factor(n):
    f = {}; d = 2
    while d*d <= n:
        while n % d == 0: f[d] = f.get(d, 0)+1; n //= d
        d += 1 if d == 2 else 2
    if n > 1: f[n] = f.get(n, 0)+1
    return f
def two_sq(p):
    for a in range(1, math.isqrt(p)+1):
        b2 = p - a*a; b = math.isqrt(b2)
        if b*b == b2 and a >= b: return (a, b)
def points(n):
    # all Gaussian integers of norm n in the sector 0<=arg<pi/2 (one per unit class)
    f = factor(n)
    pts = [(1, 0)]
    for p, e in f.items():
        a, b = two_sq(p); pi_ = (a, b); pib = (a, -b)
        new = []
        for z in pts:
            for t in range(e+1):
                w = z
                for _ in range(t): w = gmul(w, pi_)
                for _ in range(e-t): w = gmul(w, pib)
                new.append(w)
        pts = new
    out = []
    for z in pts:
        for _ in range(4):
            if z[0] > 0 and z[1] >= 0: out.append(z); break
            z = (-z[1], z[0])
    return sorted(set(out), key=lambda z: math.atan2(z[1], z[0]))
def best_cluster(n, k):
    P = points(n); d = len(P)
    ang = [math.atan2(z[1], z[0]) for z in P]
    best = None
    for i in range(d):
        j = (i + k - 1) % d
        sp = (ang[j] - ang[i]) % (math.pi/2)
        if best is None or sp < best[0]: best = (sp, i)
    sp, i = best
    cl = [P[(i+t) % d] for t in range(k)]
    return sp*n**0.25, cl
if __name__ == "__main__":
    n = int(sys.argv[1]); k = int(sys.argv[2])
    C, cl = best_cluster(n, k)
    f = factor(n)
    print("n=%d k=%d C=%.6f factors=%s" % (n, k, C, f))
    # orientation pattern: for each split prime p | n, valuation at pi=(a+bi)
    for z in cl:
        row = []
        for p, e in sorted(f.items()):
            a, b = two_sq(p); pi_ = (a, b)
            v = 0; w = z
            while True:
                q = gmul(w, (a, -b))
                if q[0] % p == 0 and q[1] % p == 0:
                    w = (q[0]//p, q[1]//p); v += 1
                else: break
            row.append(v)
        print(z, row)
