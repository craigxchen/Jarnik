"""Independent recount for Prop A (tools.md): S = Gaussian primes over split p <= P; S-integers w (all 4 units) with
R/2 < |w| <= R; squares of side isqrt(R)//2 (same cells as check_smooth2.py); T_m/T_{m-1} lower-bounds #hyperplanes.
Also an exact (5.1) check of the worst tuple in a few cells. Run: python3 -I propA_recount.py"""
import math, random
from collections import defaultdict
random.seed(3)
def primes_upto(n):
    s = bytearray([1])*(n+1); s[0:2]=b'\x00\x00'
    for i in range(2,int(n**.5)+1):
        if s[i]: s[i*i::i]=bytearray(len(s[i*i::i]))
    return [i for i in range(n+1) if s[i]]
def gauss_prime(p):
    for a in range(1, int(p**.5)+1):
        b2 = p - a*a; b = math.isqrt(b2)
        if b*b == b2: return (a, b)
def mul(x, y): return (x[0]*y[0]-x[1]*y[1], x[0]*y[1]+x[1]*y[0])
def run(P, R, m):
    X = R*R
    ps = [p for p in primes_upto(P) if p % 4 == 1]
    gps = []
    for p in ps:
        a, b = gauss_prime(p); gps += [(a, b), (a, -b)]
    gps.sort(key=lambda g: g[0]**2+g[1]**2)
    out = []
    stack = [(0, (1, 0))]
    while stack:
        k, w = stack.pop()
        n = w[0]**2 + w[1]**2
        if 4*n > X: out.append(w)
        for i in range(k, len(gps)):
            g = gps[i]
            if n*(g[0]**2+g[1]**2) > X: break
            stack.append((i, mul(w, g)))
    pts = []
    for w in out:
        x, y = w
        pts += [(x, y), (-y, x), (-x, -y), (y, -x)]
    side = math.isqrt(R)//2
    cells = defaultdict(list)
    for w in pts: cells[(w[0]//side, w[1]//side)].append(w)
    Tm = sum(math.perm(len(v), m) for v in cells.values()); Tm1 = sum(math.perm(len(v), m-1) for v in cells.values())
    # worst-case (5.1) exponent over all m-tuples in a cell, via norms and max diff in that cell (gcd ignored => upper bound on H, conservative)
    dworst = float('inf')
    for v in random.sample([v for v in cells.values() if len(v) >= m], 200):
        H = max(a*a+b*b for a, b in v); Nmin = min(a*a+b*b for a, b in v)
        D2 = max((p[0]-q[0])**2+(p[1]-q[1])**2 for p in v for q in v)
        # with g = gcd: H_x = H/N(g), N(x_i) >= Nmin/N(g), diffs/|g|: exponent is scale-free except the H normalisation; use g=1 (worst)
        lg = math.log(H) + (m-1)*math.log(D2) - m*math.log(H) - m*math.log(Nmin)
        dworst = min(dworst, -lg/math.log(H) - m)
    print(f"P={P} s={len(ps)} R={R} m={m}: U={len(pts)} T_m/T_(m-1)={Tm/Tm1:.2f}; worst-case delta_eff over 200 random cells (g=1) = {dworst:.3f}, window (1,{(m-1)/2})")
for P, R, m in [(100, 10**4, 4), (200, 10**4, 4)]:
    run(P, R, m)
