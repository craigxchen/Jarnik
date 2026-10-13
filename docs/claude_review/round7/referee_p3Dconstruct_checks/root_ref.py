# Referee's independent check of Lemma 3.2 / Prop 3.3 (root-of-unity parity repair).
# Own Z[i]/Q arithmetic; does not import the author's code.
import random, math
random.seed(2026)

def is_prime(n):
    if n < 2: return False
    i = 2
    while i*i <= n:
        if n % i == 0: return False
        i += 1
    return True

def legendre(a, q):
    a %= q
    if a == 0: return 0
    return 1 if pow(a, (q-1)//2, q) == 1 else -1

def mul(x, y, Q): return ((x[0]*y[0]-x[1]*y[1]) % Q, (x[0]*y[1]+x[1]*y[0]) % Q)
def conj(x, Q): return (x[0] % Q, (-x[1]) % Q)
def norm(x, Q): return (x[0]*x[0]+x[1]*x[1]) % Q
def gpow(x, e, Q):
    r = (1, 0)
    if e < 0:
        x = ginv(x, Q); e = -e
    while e:
        if e & 1: r = mul(r, x, Q)
        x = mul(x, x, Q); e >>= 1
    return r
def ginv(x, Q):
    n = norm(x, Q); ni = pow(n, -1, Q)
    c = conj(x, Q); return (c[0]*ni % Q, c[1]*ni % Q)

def T_elements(Q):
    return [(a, b) for a in range(Q) for b in range(Q) if (a*a+b*b) % Q == 1]

def order(x, Q, m):
    # order of x dividing m
    o = m
    for p in factor(m):
        while o % p == 0 and gpow(x, o//p, Q) == (1, 0):
            o //= p
    return o

def factor(n):
    fs = []; d = 2
    while d*d <= n:
        if n % d == 0:
            fs.append(d)
            while n % d == 0: n //= d
        d += 1
    if n > 1: fs.append(n)
    return fs

def v2(n):
    k = 0
    while n % 2 == 0: n //= 2; k += 1
    return k

def gauss_mul(x, y): return (x[0]*y[0]-x[1]*y[1], x[0]*y[1]+x[1]*y[0])
def gauss_conj(x): return (x[0], -x[1])
def gauss_pow(x, e):
    r = (1, 0)
    for _ in range(e): r = gauss_mul(r, x)
    return r

pp = []
for q in range(3, 401, 2):
    if is_prime(q):
        a = 1
        while q**a <= 400:
            pp.append((q, a)); a += 1
print("odd prime powers <= 400:", len(pp))

primes1 = [p for p in range(5, 20000) if is_prime(p) and p % 4 == 1]
# Gaussian primes of odd norm (split): pi = (x,y) with x^2+y^2 = p
gprimes = []
for p in primes1[:400]:
    for x in range(1, int(math.isqrt(p))+1):
        y2 = p - x*x; y = math.isqrt(y2)
        if y*y == y2 and y > 0:
            gprimes.append((x, y)); break

total_coll = 0; fails = 0
stats = {}
for (q, a) in pp:
    Q = q**a
    chi = 1 if q % 4 == 1 else -1
    m = (q - chi) * q**(a-1)
    T = T_elements(Q)
    assert len(T) == m, (Q, len(T), m)
    # generator
    g = None
    for t in T:
        if order(t, Q, m) == m: g = t; break
    assert g is not None
    # discrete log table
    dlog = {}
    cur = (1, 0)
    for e in range(m):
        dlog[cur] = e; cur = mul(cur, g, Q)
    k = v2(m); o = m >> k
    f = gpow(g, o, Q)               # 2-Sylow generator
    assert order(f, Q, m) == 2**k and dlog[f] % 2 == 1
    if k == 2: assert f in [(0, 1), (0, Q-1)]
    assert (k == 2) == (q % 8 in (3, 5))
    assert (k == 3) == (q % 16 in (7, 9))
    # every odd-label element has order divisible by 2^k
    for e in range(1, m, 2):
        assert (m // math.gcd(m, e)) % (2**k) == 0
    iota = dlog[(0, 1)]
    T4 = {(iota*s) % m for s in range(4)}
    # nu(z/conj z) = (Norm z / q)
    for _ in range(50):
        z = (random.randrange(Q), random.randrange(Q))
        if norm(z, Q) % q == 0: continue
        t = mul(z, ginv(conj(z, Q), Q), Q)
        assert (dlog[t] % 2 == 0) == (legendre(norm(z, Q), q) == 1)
    # random columns: zeta_j Gaussian primes with norm coprime to q; p_j random primes = 1 mod 4, != q
    for trial in range(60):
        cols = []
        while len(cols) < 5:
            zeta = random.choice(gprimes)
            if (zeta[0]**2 + zeta[1]**2) % q == 0: continue
            p = random.choice(primes1)
            if p == q: continue
            cols.append((zeta, p))
        labs = []; eps = []
        for zeta, p in cols:
            Nz = zeta[0]**2 + zeta[1]**2
            e_j = 1 if legendre(Nz, q) != legendre(p, q) else 0
            zr = (zeta[0] % Q, zeta[1] % Q)
            t = mul(mul(zr, ginv(conj(zr, Q), Q), Q), gpow(f, e_j, Q), Q)
            # parity of label equals [(p/q) = -1]
            assert dlog[t] % 2 == (1 if legendre(p, q) == -1 else 0)
            # realise r with r/conj r = t and Norm r = p mod Q (search: r = u*s, s rational)
            u = None
            for c in [(1, 0), (0, 1), (1, 1), (2, 1), (1, 2), (3, 1)]:
                tc = mul(t, conj(c, Q), Q)
                cand = ((c[0] + tc[0]) % Q, (c[1] + tc[1]) % Q)
                if norm(cand, Q) % q:
                    u = cand; break
            assert u is not None and mul(u, ginv(conj(u, Q), Q), Q) == t
            ratio = p * pow(norm(u, Q), -1, Q) % Q
            s = next(s for s in range(1, Q) if s*s % Q == ratio)
            r = (u[0]*s % Q, u[1]*s % Q)
            assert norm(r, Q) == p % Q and mul(r, ginv(conj(r, Q), Q), Q) == t
            labs.append(dlog[t]); eps.append(e_j)
        # random monomials
        for _ in range(400):
            v = [random.randint(-3, 3) for _ in range(5)]
            if all(x == 0 for x in v): continue
            L = sum(vi*li for vi, li in zip(v, labs)) % m
            # exact gamma_v = alpha conj(beta)
            alpha = (1, 0); beta = (1, 0)
            for vi, (zeta, p) in zip(v, cols):
                if vi > 0: alpha = gauss_mul(alpha, gauss_pow(zeta, vi))
                elif vi < 0: beta = gauss_mul(beta, gauss_pow(zeta, -vi))
            gam = gauss_mul(alpha, gauss_conj(beta))
            im = gauss_pow(gam, 2**k)[1]
            coll = L in T4
            E = sum(vi*ei for vi, ei in zip(v, eps))
            gm = (gam[0] % Q, gam[1] % Q)
            Zr = mul(gm, ginv(conj(gm, Q), Q), Q)
            # Prop 3.3(3) for k=2: collision iff gamma/conj gamma is a unit mod Q
            if k == 2:
                assert coll == (dlog[Zr] in T4)
            if coll:
                total_coll += 1
                if im % Q != 0: fails += 1
                if im == 0: print("Im=0 at", Q, v)
    # Prop 3.3(4): Z_v = 1 -> collides at every k=2 prime with unit i^{+-E}
    stats[Q] = k
print("collisions tested:", total_coll, "failures of Q | Im(gamma^(2^k)):", fails)
print("ALL OK" if fails == 0 else "FAIL")
