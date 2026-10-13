"""Exact check of Proposition D (Construction D): parity, and m_v | Im(Z_v^16) for all ||v||_1 <= 2,
at every prime power q^a <= X. Also: norm-3 behaviour (evidence) and comparison with random labels.
Usage: python3 check_D.py [X] [r] [seed]"""
import sys, math, random, itertools
from tq import *

X = int(sys.argv[1]) if len(sys.argv) > 1 else 400
r = int(sys.argv[2]) if len(sys.argv) > 2 else 12
seed = int(sys.argv[3]) if len(sys.argv) > 3 else 1
rng = random.Random(seed)
S0 = {17, 31, 47}
MOD0 = 4*17*31*47

def two_squares(p):
    # p = 1 mod 4 prime: return (a,b) with a^2+b^2 = p
    for c in range(2, p):
        if pow(c, (p-1)//2, p) == p-1: break
    x = pow(c, (p-1)//4, p)
    a, b = p, x
    while b*b > p: a, b = b, a % b
    c2 = p - b*b; d = math.isqrt(c2); assert d*d == c2
    return (b, d)

# profile primes p_j > X, = 1 mod 4
ps = set()
while len(ps) < r:
    p = rng.randrange(10**4, 10**6)
    if p % 4 == 1 and is_prime(p) and p > X: ps.add(p)
ps = sorted(ps)
# sub-genuine zeta_j: next prime p'_j = p_j mod MOD0, p'_j not in ps, distinct
used = set(ps); zp = []
for p in ps:
    c = p + MOD0
    while not (is_prime(c) and c not in used): c += MOD0
    used.add(c); zp.append(c)
zeta = [two_squares(c) for c in zp]

def leg(a, q):
    t = pow(a % q, (q-1)//2, q)
    return -1 if t == q-1 else t

def nu(j, q): return 1 if leg(ps[j], q) == -1 else 0

def iota_map(N, b):
    """Lemma 5.1: dict x in Z/(N/4) -> t in Z/N."""
    n4 = N//4
    targets = []
    seen = set()
    for t in range(N):
        if t % 2 != b or t in (0, N//2) or t in seen: continue
        seen.add(t); seen.add((-t) % N); targets.append(t)
    srcs = []; seen = set()
    for x in range(n4):
        if x in seen: continue
        seen.add(x); seen.add((-x) % n4); srcs.append(x)
    assert len(srcs) <= len(targets), (N, len(srcs), len(targets))
    mp = {}
    for x, t in zip(srcs, targets):
        mp[x] = t
        if (-x) % n4 != x: mp[(-x) % n4] = (-t) % N
    return mp

def crt(a1, m1, a2, m2):
    return (a1 + m1 * ((a2 - a1) * pow(m1, -1, m2) % m2)) % (m1*m2) if m2 > 1 else a1 % m1

moduli = []   # (q, A, TQ, labels list)
for q in primes_upto(X):
    if q == 2: continue
    A = 1
    while q**(A+1) <= X: A += 1
    T = TQ(q, A); m = T.m; qc = q - chi4(q)
    L = [T.dlog_ratio(z) for z in zeta]
    nus = [nu(j, q) for j in range(r)]
    if q % 8 in (3, 5):
        ell = [(L[j] + T.iota*((L[j] - nus[j]) % 2)) % m for j in range(r)]; case = 1
    elif q in S0:
        ell = L[:]; case = 2
        assert all(L[j] % 2 == nus[j] for j in range(r)), "S0 parity"
    else:
        b = (qc & -qc).bit_length() - 1
        qp = q**(A-1)
        if b == 3:
            o = qc // 8; case = 3
            ell = []
            for j in range(r):
                x1 = L[j] % qc
                p1 = crt(nus[j], 8, x1 % o, o)
                ell.append(crt(p1, qc, L[j] % qp, qp) if qp > 1 else p1)
        else:
            case = 4; Nq = qc // 4
            maps = [iota_map(Nq, 0), iota_map(Nq, 1)]
            ell = []
            for j in range(r):
                x = (L[j] % qc) % (Nq//4)
                t = maps[nus[j]][x]
                ell.append(crt(t, qc, L[j] % qp, qp) if qp > 1 else t)
    assert all(ell[j] % 2 == nus[j] for j in range(r)), ("parity", q)
    moduli.append((q, A, T, ell, L, case))

def gpow_exact(z, e):
    x, y = 1, 0
    for _ in range(e): x, y = x*z[0]-y*z[1], x*z[1]+y*z[0]
    return (x, y)

def Zv(v):
    x, y = 1, 0
    for j, c in enumerate(v):
        z = zeta[j] if c > 0 else (zeta[j][0], -zeta[j][1])
        for _ in range(abs(c)): x, y = x*z[0]-y*z[1], x*z[1]+y*z[0]
    return (x, y)

def mv(v, labs_of):
    tot = 1
    for (q, A, T, ell, L, case) in moduli:
        lab = labs_of(q, ell, L)
        s = sum(c*lab[j] for j, c in enumerate(v))
        amax = 0
        for a in range(1, A+1):
            ma = (q - chi4(q)) * q**(a-1)
            if (s % ma) % (ma//4) == 0: amax = a
            else: break
        tot *= q**amax
    return tot

def monomials(k):
    out = set()
    for supp in itertools.combinations_with_replacement(range(r), k):
        for signs in itertools.product((1, -1), repeat=k):
            v = [0]*r
            for j, sg in zip(supp, signs): v[j] += sg
            if any(v):
                # normalise sign
                first = next(c for c in v if c)
                if first < 0: v = [-c for c in v]
                if sum(map(abs, v)) == k: out.add(tuple(v))
    return sorted(out)

ours = lambda q, ell, L: ell
print(f"X={X} r={r} seed={seed}; moduli={len(moduli)}; cases:",
      {c: sum(1 for t in moduli if t[5] == c) for c in (1,2,3,4)})
fails = 0; worst = 0.0; cnt = 0
for k in (1, 2):
    for v in monomials(k):
        m = mv(v, ours); Z = Zv(v); I = gpow_exact(Z, 16)[1]
        assert I != 0
        if I % m != 0: fails += 1
        worst = max(worst, math.log(m)); cnt += 1
print(f"norm<=2: {cnt} monomials, divisibility m_v | Im(Z_v^16) failures: {fails}; max log m_v = {worst:.2f};",
      f"max_j log Norm zeta_j = {max(math.log(c) for c in zp):.2f}")
# norm 3: evidence
f3 = 0; w3 = 0.0; c3 = 0; f3case3 = 0
for v in monomials(3):
    m = mv(v, ours); Z = Zv(v); I = gpow_exact(Z, 16)[1]
    if I % m != 0: f3 += 1
    # restricted to cases 1-3 the divisibility must still hold
    m13 = 1
    for (q, A, T, ell, L, case) in moduli:
        if case == 4: continue
        s = sum(c*ell[j] for j, c in enumerate(v)); amax = 0
        for a in range(1, A+1):
            ma = (q - chi4(q)) * q**(a-1)
            if (s % ma) % (ma//4) == 0: amax = a
            else: break
        m13 *= q**amax
    if I % m13 != 0: f3case3 += 1
    w3 = max(w3, math.log(m)); c3 += 1
print(f"norm 3: {c3} monomials; m_v not dividing Im(Z_v^16): {f3} (expected > 0: the encoding ignores 3-term relations);",
      f"restricted to cases 1-3 failures: {f3case3} (must be 0); max log m_v = {w3:.2f}")
# random parity-respecting labels: max pair demand
rr = random.Random(seed+99)
rand = {}
for (q, A, T, ell, L, case) in moduli:
    m = T.m
    rand[q] = [ (2*rr.randrange(m//2) + nu(j, q)) % m for j in range(r)]
randlab = lambda q, ell, L: rand[q]
wr = max(math.log(mv(v, randlab)) for v in monomials(2))
print(f"random parity-respecting labels: max log m_v over norm 2 = {wr:.2f}")
print("ALL OK" if fails == 0 and f3case3 == 0 else "FAIL")
