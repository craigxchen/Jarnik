# Judge's own test of VIII.1 (M <= 1 + omega(N') + K_C) on integer circles.
# Exact Gaussian arithmetic; exact rational affine rank. No imports from Astra or referees.
import sys, math
from fractions import Fraction

NMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 10**6
CMAX = float(sys.argv[2]) if len(sys.argv) > 2 else 6.0
EXTRA = sys.argv[3:]  # optional explicit N values

# ---------- Gaussian integers as (a,b) ----------
def gmul(x, y): return (x[0]*y[0]-x[1]*y[1], x[0]*y[1]+x[1]*y[0])
def gnorm(x): return x[0]*x[0]+x[1]*x[1]
def gconj(x): return (x[0], -x[1])
def gdivexact(x, y):
    n = gnorm(y); t = gmul(x, gconj(y))
    if t[0] % n or t[1] % n: return None
    return (t[0]//n, t[1]//n)
def gdivround(x, y):
    n = gnorm(y); t = gmul(x, gconj(y))
    # nearest Gaussian integer
    return ((2*t[0]+n)//(2*n), (2*t[1]+n)//(2*n))
def gsub(x, y): return (x[0]-y[0], x[1]-y[1])
def ggcd(x, y):
    while y != (0, 0):
        q = gdivround(x, y); x, y = y, gsub(x, gmul(q, y))
    return x
def gval(z, pi):
    v = 0
    while True:
        q = gdivexact(z, pi)
        if q is None: return v
        z = q; v += 1

# ---------- primes ----------
def spf_sieve(n):
    s = list(range(n+1))
    for i in range(2, int(n**0.5)+1):
        if s[i] == i:
            for j in range(i*i, n+1, i):
                if s[j] == j: s[j] = i
    return s
LIM = max([NMAX] + [0])
def factor_trial(n):
    f = {}; d = 2
    while d*d <= n:
        while n % d == 0: f[d] = f.get(d, 0)+1; n //= d
        d += 1
    if n > 1: f[n] = f.get(n, 0)+1
    return f
SPF = spf_sieve(NMAX)
def factor(n):
    if n <= NMAX:
        f = {}
        while n > 1:
            p = SPF[n]; f[p] = f.get(p, 0)+1; n //= p
        return f
    return factor_trial(n)
PICACHE = {}
def pi_of(p):  # fixed Gaussian prime a+bi, a>b>0
    if p in PICACHE: return PICACHE[p]
    a = 1
    while True:
        b2 = p - a*a; b = math.isqrt(b2)
        if b*b == b2 and a > b: PICACHE[p] = (a, b); return (a, b)
        a += 1
        if a*a > p: break
    # fallback any
    for a in range(1, math.isqrt(p)+1):
        b2 = p-a*a; b = math.isqrt(b2)
        if b*b == b2: PICACHE[p] = (a, b); return (a, b)

def gpow(x, k):
    r = (1, 0)
    for _ in range(k): r = gmul(r, x)
    return r

def lattice_points(N):
    f = factor(N)
    base = [(1, 0)]
    for p, e in f.items():
        if p == 2:
            base = [gmul(b, gpow((1, 1), e)) for b in base]
        elif p % 4 == 3:
            if e % 2: return []
            base = [gmul(b, (p**(e//2), 0)) for b in base]
        else:
            pi = pi_of(p); pc = gconj(pi); new = []
            for a in range(e+1):
                t = gmul(gpow(pi, a), gpow(pc, e-a))
                new += [gmul(b, t) for b in base]
            base = new
    pts = set()
    for b in base:
        z = b
        for _ in range(4):
            pts.add(z); z = gmul(z, (0, 1))
    return sorted(pts, key=lambda z: math.atan2(z[1], z[0]) % (2*math.pi))

def affine_rank(rows):
    if len(rows) <= 1: return 0
    r0 = rows[0]
    M = [[Fraction(x - y) for x, y in zip(r, r0)] for r in rows[1:]]
    rank = 0; ncol = len(r0)
    for c in range(ncol):
        piv = None
        for i in range(rank, len(M)):
            if M[i][c] != 0: piv = i; break
        if piv is None: continue
        M[rank], M[piv] = M[piv], M[rank]
        for i in range(len(M)):
            if i != rank and M[i][c] != 0:
                fct = M[i][c]/M[rank][c]
                M[i] = [a - fct*b for a, b in zip(M[i], M[rank])]
        rank += 1
        if rank == len(M): break
    return rank

INERT = [l for l in range(3, 100000, 4) if all(l % d for d in range(2, math.isqrt(l)+1))]
def Sset(C, sharp):
    target = math.log(C) - (0.5*math.log(2) if sharp else 0.0)
    S = []; s = 0.0
    for l in INERT:
        if s > target: break
        S.append(l); s += math.log(l)/(l+1)
    return S

stats = dict(windows=0, viol_ast=0, viol_sharp=0, viol_identity=0, viol_1pi=0, viol_bound=0, viol_small=0)
maxdef = {}   # bucket of chord C' -> max (M-1-omega)
maxrdef = {}  # bucket -> max (M-1-affrank)
first_rdef = []  # (C', arcC, N, M, omega, rank)
minarc = {}

def analyse(N):
    pts = lattice_points(N)
    if len(pts) < 8: return
    n = len(pts); ang = [math.atan2(z[1], z[0]) % (2*math.pi) for z in pts]
    R4 = N**0.25
    for s in range(n):
        if ang[s] >= math.pi/2: continue
        for e in range(s+2, s+n):   # window pts[s..e], e exclusive end index: size>=3
            span = (ang[(e-1) % n] - ang[s]) % (2*math.pi)
            if span*R4 > CMAX: break
            W_ = [pts[k % n] for k in range(s, e)]
            test_window(N, W_, span*R4)

def test_window(N, Z, arcC):
    M = len(Z)
    g = Z[0]
    for z in Z[1:]: g = ggcd(g, z)
    Zp = [gdivexact(z, g) for z in Z]
    Np = N // gnorm(g)
    f = factor(Np)
    assert 2 not in f and all(p % 4 == 1 for p in f), (N, Np)
    primes = sorted(f)
    rows = [[gval(z, pi_of(p)) for p in primes] for z in Zp]
    chord = max(math.sqrt(gnorm(gsub(a, b))) for a in Zp for b in Zp)
    Cp = chord / Np**0.25
    stats['windows'] += 1
    # exact identity |g_ij|^2 = prod p^(e_p - |a_ip - a_jp|) and (1+i) | b_ij
    for i in range(M):
        for j in range(i+1, M):
            gij = ggcd(Zp[i], Zp[j])
            val = 1
            for k, p in enumerate(primes):
                val *= p**(f[p] - abs(rows[i][k]-rows[j][k]))
            if gnorm(gij) != val: stats['viol_identity'] += 1
            b = gdivexact(gsub(Zp[i], Zp[j]), gij)
            if gdivexact(b, (1, 1)) is None: stats['viol_1pi'] += 1
    om = len(primes)
    rk = affine_rank(rows)
    for sharp, key in ((False, 'viol_ast'), (True, 'viol_sharp')):
        S = Sset(Cp, sharp) if Cp > (1.0 if not sharp else math.sqrt(2)) else []
        enl = [list(r) for r in rows]
        for l in S:
            classes = sorted(set(((z[0] % l), (z[1] % l)) for z in Zp))
            idx = {c: t for t, c in enumerate(classes)}
            for zi, z in enumerate(Zp):
                oh = [0]*len(classes); oh[idx[(z[0] % l, z[1] % l)]] = 1
                enl[zi] += oh
        if affine_rank(enl) != M-1: stats[key] += 1
        if sharp and M > 1 + om + sum(S): stats['viol_bound'] += 1
    if Cp < math.sqrt(2) and rk != M-1: stats['viol_small'] += 1
    bk = round(Cp*4)/4
    maxdef[bk] = max(maxdef.get(bk, -99), M-1-om)
    maxrdef[bk] = max(maxrdef.get(bk, -99), M-1-rk)
    if N >= 10**4 and rk < M-1:
        dd = M-1-rk
        if arcC < minarc.get(dd, (99,))[0]: minarc[dd] = (round(arcC, 4), N, M, om)
    if rk < M-1: first_rdef.append((round(Cp, 4), round(arcC, 4), N, M, om, rk))

for N in range(5, NMAX+1):
    analyse(N)
for x in EXTRA:
    analyse(int(x))
print('NMAX', NMAX, 'CMAX', CMAX)
print(stats)
print('max M-1-omega by chord bucket:', sorted(maxdef.items()))
print('max M-1-affrank by chord bucket:', sorted(maxrdef.items()))
first_rdef.sort()
print('smallest-chord rank deficiencies:', first_rdef[:8])
big = [t for t in first_rdef if t[2] >= 10**4]
print('smallest-chord rank deficiencies, N>=1e4:', big[:5])
big.sort(key=lambda t: t[1]); print('smallest-arc rank deficiencies, N>=1e4:', big[:5])

print('min arc constant of affine-rank deficiency d (N>=1e4):', sorted(minarc.items()))
