# Exact test of the 2-adic block refinement of VIII.1 on integer circles.
# For each circle N <= NMAX, all windows of consecutive lattice points with arc <= ARC*sqrt(R):
# normalise by Gaussian gcd, compute chord constant C', exponent rows, residue classes mod (1+i)^k.
# Check: if C' < 2^{(1/2) sum_{j<=k} 1/f_j} then (exponent rows | one-hot class mod (1+i)^k) affinely independent.
import sys, math
from fractions import Fraction
from sympy import factorint

NMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 100000
ARC = float(sys.argv[2]) if len(sys.argv) > 2 else 3.5

def gmul(a, b): return (a[0]*b[0]-a[1]*b[1], a[0]*b[1]+a[1]*b[0])
def gnorm(a): return a[0]*a[0]+a[1]*a[1]
def gdivmod(a, b):
    n = gnorm(b); c = gmul(a, (b[0], -b[1]))
    # rounded quotient
    q = ((2*c[0] + n)//(2*n), (2*c[1] + n)//(2*n))
    r = (a[0]-gmul(q, b)[0], a[1]-gmul(q, b)[1])
    return q, r
def ggcd(a, b):
    while b != (0, 0):
        _, r = gdivmod(a, b); a, b = b, r
    return a
def gexact_div(a, b):
    q, r = gdivmod(a, b); assert r == (0, 0); return q
def divides(b, a):
    _, r = gdivmod(a, b); return r == (0, 0)

def two_sq(p):
    for a in range(1, int(math.isqrt(p)) + 1):
        b2 = p - a*a; b = math.isqrt(b2)
        if b*b == b2: return (a, b)

def red_class(z, k):
    # class of z modulo (1+i)^k, canonical
    x, y = z; h = k // 2; M = 2 ** h
    if k % 2 == 0: return (x % M, y % M)
    t = y // M; x2, y2 = x - t*M, y - t*M
    return (x2 % (2*M), y2)

def affine_rank(rows):
    base = rows[0]
    mat = [[Fraction(r[j] - base[j]) for j in range(len(base))] for r in rows[1:]]
    rank = 0; ncol = len(base)
    for c in range(ncol):
        piv = None
        for i in range(rank, len(mat)):
            if mat[i][c] != 0: piv = i; break
        if piv is None: continue
        mat[rank], mat[piv] = mat[piv], mat[rank]
        for i in range(len(mat)):
            if i != rank and mat[i][c] != 0:
                f = mat[i][c] / mat[rank][c]
                mat[i] = [mat[i][j] - f*mat[rank][j] for j in range(ncol)]
        rank += 1
    return rank

fib = {1: 1, 2: 2, 3: 4, 4: 4, 5: 8, 6: 8}
thr = {}
acc = 0
for k in range(1, 7):
    acc += Fraction(1, fib[k]); thr[k] = 2 ** (float(acc) / 2)

stats = {k: [0, 0] for k in range(1, 7)}  # tested, failures
bare_dep_below = {}
nwin = 0
for N in range(5, NMAX + 1):
    f = factorint(N)
    if any(p % 4 == 3 and e % 2 for p, e in f.items()): continue
    # lattice points
    pts = []
    for x in range(0, math.isqrt(N) + 1):
        y2 = N - x*x; y = math.isqrt(y2)
        if y*y == y2:
            for sx in (1, -1):
                for sy in (1, -1):
                    pts.append((sx*x, sy*y))
    pts = sorted(set(pts), key=lambda z: math.atan2(z[1], z[0]))
    n = len(pts)
    if n < 6: continue
    R = math.sqrt(N); L = ARC * math.sqrt(R)
    ang = [math.atan2(z[1], z[0]) for z in pts]
    for i in range(n):
        j = i
        win = [pts[i]]
        while True:
            jj = (j + 1) % n
            if jj == i: break
            d = (ang[jj] - ang[i]) % (2*math.pi)
            if R * d > L: break
            win.append(pts[jj]); j = jj
        if len(win) < 3: continue
        # all sub-windows starting at i (consecutive prefixes) of size >= 3
        for M in range(3, len(win) + 1):
            W_ = win[:M]
            g = W_[0]
            for z in W_[1:]: g = ggcd(g, z)
            zs = [gexact_div(z, g) for z in W_]
            Np = gnorm(zs[0])
            Cp = max(math.sqrt(gnorm((a[0]-b[0], a[1]-b[1]))) for a in zs for b in zs) / Np ** 0.25
            fp = factorint(Np)
            assert all(p % 4 == 1 for p in fp), (N, fp)
            rows = []
            for z in zs:
                r = []
                for p in sorted(fp):
                    pi = two_sq(p); pi = (pi[0], pi[1]); v = 0; w = z
                    while divides(pi, w): w = gexact_div(w, pi); v += 1
                    r.append(v)
                rows.append(r)
            nwin += 1
            s = len(fp)
            if Cp < 2 and affine_rank(rows) < M - 1:
                bare_dep_below[(N, M)] = Cp
            for k in range(1, 7):
                if Cp < thr[k]:
                    classes = sorted(set(red_class(z, k) for z in zs))
                    assert len(classes) <= fib[k]
                    ext = [r + [1 if red_class(z, k) == c else 0 for c in classes] for r, z in zip(rows, zs)]
                    stats[k][0] += 1
                    if affine_rank(ext) < M - 1:
                        stats[k][1] += 1
                        print("FAIL", k, N, M, Cp)
                    if M > s + fib[k]:
                        print("SIZEFAIL", k, N, M, s)
print("windows", nwin)
for k in range(1, 7):
    print("k", k, "threshold C' <", round(thr[k], 5), "cost", fib[k]-1, "tested", stats[k][0], "failures", stats[k][1])
print("bare-row dependences at C'<2:", len(bare_dep_below))
