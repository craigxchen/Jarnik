"""Checks for lit18/tools.md (evidence + exact identities; scratch only).

(A) Lines in X_M through every point of the diagonal line Y (so the Seshadri constant eps_Y(H) <= 1).
(B) Pigeonhole pairs of Gaussian primes whose ratio rho = eta*conj(u)/u has gcd-ratio
    -log|rho-1| / h(rho) > 1/2 (h(rho) = (1/2) log(p_i p_j)): these force superlinear exceptional
    sets in any S-uniform two-function gcd theorem at epsilon < 1/2.
(C) Smooth-number tuples: m S-integers (S = Gaussian primes over split p <= P) in a common square of
    side sqrt(R)/2 solve the round-8 system (5.1) with delta close to (m-1)/2 (dense-S regime).
Run: python3 -I check_tools.py
"""
import math, itertools, sys
from fractions import Fraction

# ---------- (A) ----------
def check_lines(M=6):
    import sympy as sp
    u, v, c, t = sp.symbols('u v c t', nonzero=True)
    J1 = range(0, 2); J2 = range(2, M)
    a = [0 if j in J1 else c / v for j in range(M)]
    b = [c / u if j in J1 else 0 for j in range(M)]
    vals = [sp.expand((u + t * a[j]) * (v + t * b[j])) for j in range(M)]
    ok = all(sp.simplify(vals[j] - vals[0]) == 0 for j in range(M))
    # the line leaves Y for t != 0: u_0 + t a_0 = u but u_2 + t a_2 = u + t c/v
    leaves = sp.simplify((u + t * a[2]) - (u + t * a[0])) != 0
    print(f"(A) M={M}: line p + t(a;b) lies on X_M: {ok}; not contained in Y: {leaves}")

# ---------- Gaussian helpers ----------
def primes_upto(n):
    s = bytearray([1]) * (n + 1); s[0:2] = b'\x00\x00'
    for i in range(2, int(n ** 0.5) + 1):
        if s[i]:
            s[i * i::i] = bytearray(len(s[i * i::i]))
    return [i for i in range(n + 1) if s[i]]

def gauss_prime(p):  # p = 1 mod 4 -> (a,b), a>b>0, a^2+b^2=p
    for b in range(1, int(math.isqrt(p)) + 1):
        a2 = p - b * b; a = math.isqrt(a2)
        if a * a == a2 and a > b:
            return (a, b)
    raise ValueError

def gmul(x, y):
    return (x[0] * y[0] - x[1] * y[1], x[0] * y[1] + x[1] * y[0])

def gconj(x):
    return (x[0], -x[1])

def gnorm(x):
    return x[0] * x[0] + x[1] * x[1]

def gdivmod(x, y):
    n = gnorm(y); num = gmul(x, gconj(y))
    q = ((2 * num[0] + n) // (2 * n), (2 * num[1] + n) // (2 * n))
    r = (x[0] - gmul(q, y)[0], x[1] - gmul(q, y)[1])
    return q, r

def ggcd(x, y):
    while y != (0, 0):
        _, r = gdivmod(x, y); x, y = y, r
    return x

def gexact_div(x, y):
    q, r = gdivmod(x, y); assert r == (0, 0); return q

# ---------- (B) ----------
def check_pairs(slist=(50, 100, 200, 400, 800, 1600), thr=0.5):
    ps = [p for p in primes_upto(200000) if p % 4 == 1]
    units = [(1, 0), (0, 1), (-1, 0), (0, -1)]
    print(f"(B) pairs (i<j among first s split primes) with gcd-ratio -log|rho-1|/h(rho) > {thr}")
    for s in slist:
        P = ps[:s]; G = [gauss_prime(p) for p in P]
        th = [math.atan2(g[1], g[0]) for g in G]
        cnt = 0; best = 0.0
        for i in range(s):
            for j in range(i + 1, s):
                # rho = eta * conj(u)/u, u = pi_i conj(pi_j); arg(rho) = -2(th_i - th_j) + arg(eta)
                ang = -2.0 * (th[i] - th[j])
                d = min(abs((ang - k * math.pi / 2 + math.pi) % (2 * math.pi) - math.pi) for k in range(4))
                dist = 2 * math.sin(d / 2)  # |rho - 1|
                h = 0.5 * math.log(P[i] * P[j])
                r = -math.log(dist) / h if dist > 0 else float('inf')
                best = max(best, r)
                if r > thr:
                    cnt += 1
        print(f"    s={s:5d}  P_s={P[-1]:7d}  count={cnt:7d}  count/s={cnt/s:7.3f}  max ratio={best:.3f}")

# ---------- (C) ----------
def s_integers(P, X):
    """All Gaussian integers (up to the 4 units) with every prime factor over a split p <= P, norm <= X."""
    gp = []
    for p in primes_upto(P):
        if p % 4 == 1:
            a, b = gauss_prime(p); gp += [((a, b), p), ((a, -b), p)]
    out = []
    def rec(k, cur, n):
        out.append(cur)
        for idx in range(k, len(gp)):
            g, p = gp[idx]
            if n * p > X:
                continue
            rec(idx, gmul(cur, g), n * p)
    rec(0, (1, 0), 1)
    units = [(1, 0), (0, 1), (-1, 0), (0, -1)]
    allw = [gmul(w, e) for w in out for e in units]
    return allw, sum(1 for p in primes_upto(P) if p % 4 == 1)

def check_smooth(P=200, R=10 ** 4, m=4, max_report=5):
    X = R * R
    W, s = s_integers(P, X)
    side = max(1, int(math.isqrt(R) // 2))
    cells = {}
    for w in W:
        if gnorm(w) * 4 < X:  # keep the annulus R/2 < |w| <= R
            continue
        cells.setdefault((w[0] // side, w[1] // side), []).append(w)
    U = sum(len(v) for v in cells.values())
    tm = 0; tm1 = 0; dmin = float('inf'); dmax = 0.0; shown = 0
    proj = set()
    for key, pts in cells.items():
        k = len(pts)
        tm1 += math.perm(k, m - 1)
        if k < m:
            continue
        for tup in itertools.permutations(pts, m):
            g = tup[0]
            for w in tup[1:]:
                g = ggcd(g, w)
            x = [gexact_div(w, g) for w in tup]
            H = max(gnorm(xi) for xi in x)
            fin = -sum(math.log(gnorm(xi)) for xi in x)
            arch = math.log(gnorm(x[0])) + sum(math.log(gnorm((xi[0] - x[0][0], xi[1] - x[0][1]))) for xi in x[1:]) - m * math.log(H)
            delta = -(fin + arch) / math.log(H) - m
            dmin = min(dmin, delta); dmax = max(dmax, delta); tm += 1
            # projective point via normalised ratios
            proj.add(tuple((Fraction(gmul(xi, gconj(x[0]))[0], gnorm(x[0])), Fraction(gmul(xi, gconj(x[0]))[1], gnorm(x[0]))) for xi in x[1:]))
            if shown < max_report:
                print(f"    tuple {tup}  H={H}  delta_eff={delta:.3f}")
                shown += 1
    print(f"(C) P={P} (s={s} split primes), R={R}, m={m}: S-integers in annulus U={U}, cells={len(cells)}, "
          f"ordered m-tuples in common cells={tm}, distinct projective points={len(proj)}, "
          f"ordered (m-1)-tuples={tm1}, delta_eff range=[{dmin:.3f},{dmax:.3f}] (need > 1 and < {(m-1)/2})")
    if tm1:
        print(f"    hyperplane lower bound #proj/#(m-1)-tuples = {len(proj)/tm1:.4f}")

if __name__ == '__main__':
    check_lines(6)
    check_pairs()
    for (P, R, m) in [(100, 10 ** 4, 3), (200, 10 ** 4, 4), (400, 3 * 10 ** 4, 4)]:
        check_smooth(P, R, m)
