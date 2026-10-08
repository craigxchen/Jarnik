"""Referee check R3/R4 on actual circles (own Gaussian-integer code, independent of gauss.py).

R3 (Lemma 4.1, alpha=1/2): for every subset S (2<=|S|<=7) of lattice points lying on one arc of
   length <= C sqrt(R) of x^2+y^2=N (=R^2), with d = gcd of all coordinates, R'=R/d,
   mx = max|coord|/d, md = max_j max|coord(z_j-z_1)|/d, ch2 = max squared chord (original):
     (a) d^2 <= ch2                                   [d <= |z_i - z_j|]
     (c) 2 mx^2 ch2 d >= md^2 N                       [lambda_Y >= (1/2)log R' - log C_eff - (1/2)log 2,
                                                       C_eff = sqrt(ch2)/sqrt(R)]  (exact integers)
   and the consequence used in Thm 5.1 (M=7): Phi = 5 lambda_Y - 2 h >= (1/2) log R' - 5 log C_eff - (5/2) log 2.
R4 (Prop 8.1 identity): for all pairs z != w on the circle: g = gcd(z,w), u = z/g, w = g*eta*conj(u),
   (s,t) per eta, s^2+t^2 in {n,2n}, gcd(s,t)=1, s != 0, |z-w|^2 >= 2|g|^2 s^2.
Also Machin's pair on 13^8.
"""
import random, math
from math import gcd, log
from itertools import combinations

def gmul(a, b): return (a[0]*b[0] - a[1]*b[1], a[0]*b[1] + a[1]*b[0])
def gconj(a): return (a[0], -a[1])
def gnorm(a): return a[0]*a[0] + a[1]*a[1]
def gdivmod(a, b):
    n = gnorm(b); num = gmul(a, gconj(b))
    q = ((2*num[0] + n) // (2*n), (2*num[1] + n) // (2*n))   # rounded quotient
    qb = gmul(q, b); r = (a[0]-qb[0], a[1]-qb[1])
    return q, r
def ggcd(a, b):
    while b != (0, 0):
        _, r = gdivmod(a, b); a, b = b, r
    return a
def gexact(a, b):
    q, r = gdivmod(a, b); assert r == (0, 0); return q

def prime_pi(p):  # p = 1 mod 4: return Gaussian prime of norm p
    for x in range(1, math.isqrt(p) + 1):
        y2 = p - x*x; y = math.isqrt(y2)
        if y*y == y2: return (x, y)
    raise ValueError

def circle_points(fact):
    """fact: dict prime->exponent for N. returns all Gaussian z with |z|^2 = N (empty if impossible)."""
    base = [(1, 0)]
    for p, e in fact.items():
        if p == 2:
            t = (1, 0)
            for _ in range(e): t = gmul(t, (1, 1))
            base = [gmul(b, t) for b in base]
        elif p % 4 == 3:
            if e % 2: return []
            base = [(b[0]*p**(e//2), b[1]*p**(e//2)) for b in base]
        else:
            pi = prime_pi(p); pib = gconj(pi); new = []
            for a in range(e + 1):
                t = (1, 0)
                for _ in range(a): t = gmul(t, pi)
                for _ in range(e - a): t = gmul(t, pib)
                new += [gmul(b, t) for b in base]
            base = new
    pts = set()
    for b in base:
        z = b
        for _ in range(4):
            pts.add(z); z = gmul(z, (0, 1))
    return sorted(pts, key=lambda z: math.atan2(z[1], z[0]))

def check_R3(pts, N, C, stats, maxsub=7):
    R = math.sqrt(N); L = C * math.sqrt(R); n = len(pts)
    ang = [math.atan2(z[1], z[0]) for z in pts]
    # windows of consecutive points (cyclic) with angular span*R <= L
    for i in range(n):
        win = [pts[i]]
        for k in range(1, n):
            j = (i + k) % n
            span = (ang[j] - ang[i]) % (2*math.pi)
            if span * R <= L: win.append(pts[j])
            else: break
        if len(win) < 2: continue
        for r in range(2, min(maxsub, len(win)) + 1):
            for S in combinations(win, r):
                if S[0] != pts[i]: continue   # each subset once (anchored at first point)
                coords = [c for z in S for c in z]
                d = 0
                for c in coords: d = gcd(d, c)
                Sr = [(z[0]//d, z[1]//d) for z in S]
                mx = max(abs(c) for z in Sr for c in z)
                md = max(max(abs(z[0]-Sr[0][0]), abs(z[1]-Sr[0][1])) for z in Sr)
                ch2 = max(gnorm((a[0]-b[0], a[1]-b[1])) for a, b in combinations(S, 2))
                stats['tuples'] += 1
                if not d*d <= ch2: stats['fail_a'] += 1
                if not 2*mx*mx*ch2*d >= md*md*N: stats['fail_c'] += 1
                if r == 7:
                    Rp = R / d; Ceff = math.sqrt(ch2) / math.sqrt(R)
                    lam = log(mx) - log(md); h = log(mx)
                    Phi = 5*lam - 2*h
                    bound = 0.5*log(Rp) - 5*log(Ceff) - 2.5*log(2)
                    stats['seven'] += 1
                    if Phi < bound - 1e-9: stats['fail_phi'] += 1
                    stats['min_phi_slack'] = min(stats['min_phi_slack'], Phi - bound)

def check_R4(pts, N, stats):
    for z, w in combinations(pts, 2):
        g = ggcd(z, w); u = gexact(z, g); up = gexact(w, g)
        ub = gconj(u); eta = None
        for e in [(1, 0), (-1, 0), (0, 1), (0, -1)]:
            if gmul(e, ub) == up: eta = e
        assert eta is not None
        A, B = u
        s, t = {(1, 0): (B, A), (-1, 0): (A, B), (0, 1): (A - B, A + B), (0, -1): (A + B, A - B)}[eta]
        n = gnorm(u); c0 = s*s + t*t
        ok = (c0 in (n, 2*n)) and gcd(s, t) == 1 and s != 0 and gcd(A, B) == 1
        dz2 = gnorm((z[0]-w[0], z[1]-w[1]))
        ok = ok and dz2 >= 2*gnorm(g)*s*s and N % n == 0
        stats['pairs'] += 1
        if not ok: stats['fail_pair'] += 1

if __name__ == '__main__':
    rng = random.Random(20261008)
    split = [5, 13, 17, 29, 37, 41, 53, 61, 73, 89, 97, 101, 109, 113]
    st3 = dict(tuples=0, fail_a=0, fail_c=0, seven=0, fail_phi=0, min_phi_slack=1e9)
    st4 = dict(pairs=0, fail_pair=0)
    ncirc = 0
    for trial in range(400):
        fact = {}
        for _ in range(rng.randint(2, 6)):
            p = rng.choice(split); fact[p] = fact.get(p, 0) + rng.randint(1, 3)
        if rng.random() < 0.3: fact[2] = rng.randint(1, 3)
        if rng.random() < 0.3: fact[rng.choice([3, 7, 11])] = 2
        N = 1
        for p, e in fact.items(): N *= p**e
        if N > 10**16: continue
        pts = circle_points(fact)
        assert all(gnorm(z) == N for z in pts)
        ncirc += 1
        for C in (1, 2, 4, 8):
            check_R3(pts, N, C, st3)
        if len(pts) <= 400: check_R4(pts, N, st4)
    print(f"R3: {ncirc} circles, {st3['tuples']} short-arc subtuples (C in 1,2,4,8; sizes 2..7): "
          f"fail(a)={st3['fail_a']} fail(c)={st3['fail_c']}; 7-subtuples={st3['seven']} fail(Phi bound)={st3['fail_phi']}"
          f" min slack={st3['min_phi_slack']:.3f}")
    print(f"R4: {st4['pairs']} pairs: identity failures = {st4['fail_pair']}")
    # Machin
    z = (13**4, 0); w = (28560, -239)
    print("Machin: 239^2+1 == 2*13^4:", 239**2 + 1 == 2*13**4, "; both on 13^8:", gnorm(z) == 13**8 == gnorm(w),
          "; |z-w|/sqrt(R) =", math.sqrt(gnorm((z[0]-w[0], z[1]-w[1]))) / 13**2)
    g = ggcd(z, w); print("   gcd norm =", gnorm(g), " n = R^2/|g|^2 =", 13**8 // gnorm(g), "= 13^", round(math.log(13**8 // gnorm(g), 13)))
