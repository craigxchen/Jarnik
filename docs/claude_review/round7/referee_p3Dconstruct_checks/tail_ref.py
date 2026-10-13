# Referee's independent rerun of Evidence 1 (single-column demands) - evidence only.
# COH: zeta_j distinct auxiliary Gaussian primes of norm > X, root repair (3.1) at every q.
# RAND: labels uniform in parity class. Genuine parities of random primes p_j ~ M^4.
import random, math, sys
random.seed(int(sys.argv[2]) if len(sys.argv) > 2 else 5)
def is_prime(n):
    if n < 2: return False
    if n % 2 == 0: return n == 2
    for p in [3,5,7,11,13,17,19,23,29,31,37]:
        if n % p == 0: return n == p
    d = n-1; s = 0
    while d % 2 == 0: d //= 2; s += 1
    for a in [2,3,5,7,11,13,17,19,23,29,31,37]:
        x = pow(a, d, n)
        if x in (1, n-1): continue
        for _ in range(s-1):
            x = x*x % n
            if x == n-1: break
        else: return False
    return True
def leg(a, q):
    a %= q
    return 0 if a == 0 else (1 if pow(a, (q-1)//2, q) == 1 else -1)
def mul(x, y, Q): return ((x[0]*y[0]-x[1]*y[1]) % Q, (x[0]*y[1]+x[1]*y[0]) % Q)
def gpow(x, e, Q):
    r = (1, 0)
    while e:
        if e & 1: r = mul(r, x, Q)
        x = mul(x, x, Q); e >>= 1
    return r
def ratio(z, Q):  # z / conj z mod Q
    n = (z[0]*z[0]+z[1]*z[1]) % Q; ni = pow(n, -1, Q)
    zz = mul(z, z, Q)  # z/conj z = z^2 / Norm z
    return (zz[0]*ni % Q, zz[1]*ni % Q)
M = int(sys.argv[1]); X = 3*M; r = 8*M
qs = [q for q in range(3, X+1, 2) if is_prime(q)]
info = []
for q in qs:
    A = 1
    while q**(A+1) <= X: A += 1
    Q = q**A
    chi = 1 if q % 4 == 1 else -1
    m = (q - chi)*q**(A-1)
    k = (m & -m).bit_length()-1; o = m >> k
    # 2-Sylow generator: (random z/conj z)^o of order exactly 2^k
    while True:
        z = (random.randrange(Q), random.randrange(Q))
        if (z[0]**2+z[1]**2) % q == 0: continue
        f = gpow(ratio(z, Q), o, Q)
        if gpow(f, 2**(k-1), Q) != (1, 0): break
    info.append((q, A, Q, m, k, f))
# auxiliary Gaussian primes of norm > X, distinct
aux = []; n0 = X+1
while len(aux) < r:
    if n0 % 4 == 1 and is_prime(n0):
        a = 1
        while True:
            b2 = n0 - a*a; b = math.isqrt(b2)
            if b*b == b2: aux.append((a, b)); break
            a += 1
    n0 += 1
ps = []
while len(ps) < r:
    p = random.randrange(M**4//2, M**4)
    if p % 4 == 1 and is_prime(p): ps.append(p)
coh = []; rnd = []
logM = math.log(M)
for j in range(r):
    zeta = aux[j]; Nz = zeta[0]**2+zeta[1]**2; p = ps[j]
    lc = 0.0; lr = 0.0
    for (q, A, Q, m, k, f) in info:
        e = 1 if leg(Nz, q) != leg(p, q) else 0
        t = ratio(zeta, Q)
        if e: t = mul(t, f, Q)
        t4 = gpow(t, 4, Q)
        a = 0
        for aa in range(A, 0, -1):
            Qa = q**aa
            if t4[0] % Qa == 1 % Qa and t4[1] % Qa == 0: a = aa; break
        lc += a*math.log(q)
        # RAND: uniform label of parity nu = [(p/q) = -1] mod m
        nu = 1 if leg(p, q) == -1 else 0
        lab = 2*random.randrange(m//2) + nu
        a = 0
        for aa in range(A, 0, -1):
            mq = (m // q**(A-aa))
            if lab % (mq//4) == 0: a = aa; break
        lr += a*math.log(q)
    coh.append(lc/logM); rnd.append(lr/logM)
print(f"M={M} r={r} X={X}: RAND max/mean {max(rnd):.2f}/{sum(rnd)/r:.2f}  COH max/mean {max(coh):.2f}/{sum(coh)/r:.2f}  logM/loglogM={logM/math.log(logM):.2f}")
