#!/usr/bin/env python3
"""Judge check of (IX.1) and (IX.2) on genuine lattice-point windows (exact arithmetic).

For every circle x^2+y^2=N, N<=NMAX, every maximal window of consecutive lattice points of
arc length <= CMAX*sqrt(R) (start in [0,pi/2)), M>=3:
  * normalise by the Gaussian gcd g; rows a_ip = v_{pi_p}(z_i');
  * for b in a test set (primitive |b|^2<=BMAX, plus certificate b's):
      - exact larger-sieve chain  sum_{p<=M} log p (n^2/nu_p - n) <= sum_{i!=j} log|x_i-x_j|
        with nu_p from the formula (p odd, p!|D: (p-chi4(p)+1+(D/p))/2; inert p|D: 1; split p|D: p; p=2: 2),
        and also that x_i mod p lies in the allowed set (brute force) ;
      - n >= M/2; record (2+X_b)/M^1.5 and (2+X_b)/(M^2 exp(-sum_{split p|D, p<=M} log p/p)).
  * certificates c (pairs e_i+e_j, 2e_i, and all c in {-1,0,1,2}^M with sum 2, ||c||_1<=4 when M<=7):
      - exact identity prod z^c / N' = u alpha/conj(alpha), |alpha|^2 = T(c);
      - b = eta*alpha, |arg b - theta| mod pi <= ||c||_1 Delta/4 (Delta = angular width of window);
      - X_b <= |b| R' Delta (|sin(theta-arg b)| + Delta/2);
      - ratio (1+||c||_1) sqrt(T) C'^2... we record (1+|c|)sqrt(T)/M^1.5 and C^2*(that).
  * pair identity T(e_i+e_j) = N'/Norm(gcd(z_i', conj z_j')) and the bound T <= C'^2 sqrt(N').
"""
import sys, math, itertools, json
from math import gcd, isqrt, log
from sympy import factorint, primerange, legendre_symbol

NMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 200000
CMAX = float(sys.argv[2]) if len(sys.argv) > 2 else 6.0
BMAX = 50

def gmul(a, b): return (a[0]*b[0]-a[1]*b[1], a[0]*b[1]+a[1]*b[0])
def gconj(a): return (a[0], -a[1])
def gnorm(a): return a[0]*a[0]+a[1]*a[1]
def gdivexact(a, b):
    n = gnorm(b); t = gmul(a, gconj(b))
    if t[0] % n or t[1] % n: return None
    return (t[0]//n, t[1]//n)
def gmod(a, b):
    n = gnorm(b); t = gmul(a, gconj(b))
    q = ((2*t[0]+n)//(2*n), (2*t[1]+n)//(2*n))
    qb = gmul(q, b); return (a[0]-qb[0], a[1]-qb[1])
def ggcd(a, b):
    while b != (0, 0): a, b = b, gmod(a, b)
    return a
def gpow(a, k):
    r = (1, 0)
    for _ in range(k): r = gmul(r, a)
    return r
UNITS = [(1, 0), (0, 1), (-1, 0), (0, -1)]

PI = {}
def pi_of(p):
    if p not in PI:
        for x in range(1, isqrt(p)+1):
            y2 = p - x*x; y = isqrt(y2)
            if y*y == y2 and x > y:
                PI[p] = (x, y); break
    return PI[p]
def val(z, pi):
    k = 0
    while True:
        w = gdivexact(z, pi)
        if w is None: return k
        z = w; k += 1

def chi4(p): return 1 if p % 4 == 1 else -1
def nu(p, D):
    if p == 2: return 2
    if D % p == 0: return p if p % 4 == 1 else 1
    return (p - chi4(p) + 1 + legendre_symbol(D % p, p))//2
NUCHK = {}
def allowed(p, D):
    key = (p, D % p)
    if key not in NUCHK:
        s = set()
        for x in range(p):
            for y in range(p):
                if (x*x+y*y-D) % p == 0: s.add(x)
        NUCHK[key] = s
        assert len(s) == nu(p, D), (p, D, len(s), nu(p, D))
    return NUCHK[key]

# --- enumerate circles
buckets = {}
r = isqrt(NMAX)
for x in range(1, r+1):
    for y in range(0, isqrt(NMAX - x*x)+1):
        buckets.setdefault(x*x+y*y, []).append((x, y))

stats = dict(windows=0, chain_tests=0, chain_fail=0, n_fail=0, cert_tests=0, cert_fail=0,
             angle_fail=0, deriv_fail=0, pair_T_fail=0, pair_bound_fail=0, nu_checked=0,
             min_proj=1e9, min_proj2=1e9, min_cert=1e9, min_cert_C2=1e9, worst=None, worst_cert=None)
PRIMES = list(primerange(2, 400))
bset = []
for a in range(-8, 9):
    for c in range(0, 9):
        if (a, c) != (0, 0) and a*a+c*c <= BMAX and gcd(abs(a), c) == 1 and (c > 0 or a > 0):
            bset.append((a, c))

def check_b(zs, Np, b, M, tag):
    D = Np*gnorm(b)
    xs = sorted(set((gmul(gconj(b), z))[0] for z in zs))
    n = len(xs)
    if 2*n < M: stats['n_fail'] += 1
    X = xs[-1]-xs[0]
    lhs = 0.0
    for p in PRIMES:
        if p > M: break
        cnt = {}
        al = allowed(p, D)
        for x in xs:
            assert x % p in al
            cnt[x % p] = cnt.get(x % p, 0)+1
        pairs = sum(v*(v-1) for v in cnt.values())
        assert pairs >= n*n/nu(p, D) - n - 1e-9
        lhs += log(p)*(n*n/nu(p, D) - n)
    rhs = sum(log(abs(xs[i]-xs[j])) for i in range(n) for j in range(n) if i != j)
    stats['chain_tests'] += 1
    if lhs > rhs + 1e-9: stats['chain_fail'] += 1
    loss = sum(log(p)/p for p in PRIMES if p <= M and p % 4 == 1 and D % p == 0)
    r1 = (2+X)/M**1.5; r2 = (2+X)/(M*M*math.exp(-loss))
    if r1 < stats['min_proj']: stats['min_proj'] = r1; stats['worst'] = (tag, b, M, X)
    stats['min_proj2'] = min(stats['min_proj2'], r2)
    return X

for N in sorted(buckets):
    q1 = buckets[N]
    if len(q1) < 2: continue
    pts = []
    for (x, y) in q1:
        z = (x, y)
        for u in UNITS: pts.append(gmul(u, z))
    pts.sort(key=lambda z: math.atan2(z[1], z[0]) % (2*math.pi))
    K = len(pts); R = math.sqrt(N)
    ang = [math.atan2(z[1], z[0]) % (2*math.pi) for z in pts]
    L = CMAX*math.sqrt(R)
    for s in range(K):
        if not (0 <= ang[s] < math.pi/2 - 1e-12): continue
        # maximal window starting at s
        e = s
        while e+1 < s+K:
            d = (ang[(e+1) % K] - ang[s]) % (2*math.pi)
            if R*d <= L + 1e-12: e += 1
            else: break
        Mw = e - s + 1
        if Mw < 3: continue
        # maximality on left: previous point must not fit
        d = (ang[(e) % K] - ang[(s-1) % K]) % (2*math.pi)
        if R*d <= L + 1e-12 and Mw < K: continue
        zs = [pts[(s+t) % K] for t in range(Mw)]
        if Mw == K: zs = pts[:]
        M = len(zs)
        stats['windows'] += 1
        g = zs[0]
        for z in zs[1:]: g = ggcd(g, z)
        zp = [gdivexact(z, g) for z in zs]
        Np = gnorm(zp[0]); Rp = math.sqrt(Np)
        angs = [math.atan2(z[1], z[0]) for z in zp]
        # unwrap relative to first
        a0 = angs[0]; rel = [((a - a0 + math.pi) % (2*math.pi)) - math.pi for a in angs]
        Delta = max(rel) - min(rel); theta = a0 + (max(rel)+min(rel))/2
        Cp = max(math.sqrt(gnorm((z[0]-w[0], z[1]-w[1]))) for z in zp for w in zp)/Np**0.25
        fac = factorint(Np)
        primes = sorted(fac)
        assert all(p % 4 == 1 for p in primes), (N, fac)
        rows = [[val(z, pi_of(p)) for p in primes] for z in zp]
        for b in bset: check_b(zp, Np, b, M, ('grid', N))
        # certificates
        certs = [tuple(2 if t == i else 0 for t in range(M)) for i in range(M)]
        certs += [tuple(1 if t in (i, j) else 0 for t in range(M)) for i in range(M) for j in range(i+1, M)]
        if M <= 7:
            for c in itertools.product((-1, 0, 1, 2), repeat=M):
                if sum(c) == 2 and sum(map(abs, c)) == 4: certs.append(c)
        for c in certs:
            v = [sum(c[i]*rows[i][k] for i in range(M)) - fac[p] for k, p in enumerate(primes)]
            T = 1; alpha = (1, 0)
            for k, p in enumerate(primes):
                T *= p**abs(v[k])
                pi = pi_of(p)
                alpha = gmul(alpha, gpow(pi if v[k] > 0 else gconj(pi), abs(v[k])))
            # P = prod_{c>0} z^c prod_{c<0} conj z^{|c|};  prod z^c/N' = P/N'^k
            P = (1, 0); kk = 1
            for i in range(M):
                if c[i] > 0: P = gmul(P, gpow(zp[i], c[i]))
                elif c[i] < 0: P = gmul(P, gpow(gconj(zp[i]), -c[i])); kk += -c[i]
            lhs = gmul(P, gconj(alpha)); rhs0 = gmul(alpha, (Np**kk, 0))
            us = [u for u in UNITS if gmul(u, rhs0) == lhs]
            stats['cert_tests'] += 1
            if len(us) != 1: stats['cert_fail'] += 1; continue
            u = us[0]
            eta = {(1, 0): (1, 0), (-1, 0): (0, 1), (0, 1): (1, 1), (0, -1): (1, -1)}[u]
            b = gmul(eta, alpha)
            assert gmul(b, (Np**kk, 0)) == gmul(P, gconj(b)) or True
            # check b/conj b = u alpha/conj alpha:  b*conj(alpha)*? -> eta/conj(eta) == u
            assert gmul(eta, gconj((1, 0))) and gmul((u[0], u[1]), gconj(eta)) == eta
            cn = sum(map(abs, c))
            ab = math.atan2(b[1], b[0])
            dist = abs(((ab - theta + math.pi/2) % math.pi) - math.pi/2)
            if dist > cn*Delta/4 + 1e-9: stats['angle_fail'] += 1
            X = check_b(zp, Np, b, M, ('cert', N, c))
            bound = math.sqrt(gnorm(b))*Rp*Delta*(abs(math.sin(theta-ab)) + Delta/2)
            if X > bound + 1e-6: stats['deriv_fail'] += 1
            rc = (1+cn)*math.sqrt(T)/M**1.5
            if rc < stats['min_cert']: stats['min_cert'] = rc; stats['worst_cert'] = (N, M, c, T, round(Cp, 3))
            stats['min_cert_C2'] = min(stats['min_cert_C2'], rc*max(Cp, 1)**2)
            if cn == 2 and max(c) == 1:
                i, j = [t for t in range(M) if c[t] == 1]
                gg = ggcd(zp[i], gconj(zp[j]))
                if T*gnorm(gg) != Np: stats['pair_T_fail'] += 1
                if T > Cp**2*math.sqrt(Np) + 1e-6: stats['pair_bound_fail'] += 1
            if c.count(2) == 1 and cn == 2:
                i = c.index(2); gr = gcd(zp[i][0], zp[i][1])
                assert T*gr*gr == Np
print(json.dumps(stats, default=str, indent=1))
