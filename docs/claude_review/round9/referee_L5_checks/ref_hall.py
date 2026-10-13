#!/usr/bin/env python3
"""Referee checks around L5.md Lemma 3.4 and Hypothesis IFL.

Part A (model, Lemma 3.4): Monte Carlo of Hall failure for t uniform candidates per row in an n-set,
  compared with the lemma's 1.5/n and with the exact union-bound sum; and the SUM over the moduli of
  Q_large(c) for X = 3 M^A, which Theorem 3.2 needs (Hall at every modulus simultaneously).
Part B (evidence, IFL for GENUINE labels): a random twin profile with genuine split primes p_j <= M^2,
  genuine labels ell_j(q) = log(pi_j / conj pi_j) in T_q at level 1, and S_y = {0, e_(j_1), ..., e_(j_(t-1))}
  (free coordinates).  Part (a) of IFL holds for genuine labels by the height bound (L5.md Remark 1.3);
  here we test part (b) (Hall) at every q in Q_large(c) up to a bound, and the forced Legendre parity.
"""
import random, math, sys
random.seed(11)

def is_prime(n):
    if n < 2: return False
    i = 2
    while i * i <= n:
        if n % i == 0: return False
        i += 1
    return True

def hall_ok(cands_by_row, nvals=None):
    """cands_by_row: list of sets of values; True iff a system of distinct representatives exists."""
    match = {}
    def try_row(u, seen):
        for v in cands_by_row[u]:
            if v in seen: continue
            seen.add(v)
            if v not in match or try_row(match[v], seen):
                match[v] = u
                return True
        return False
    for u in range(len(cands_by_row)):
        if not try_row(u, set()): return False
    return True

# ---------------- Part A -------------------------------------------------------------------------
def union_bound(M, n, t):
    tot = 0.0
    for s in range(2, M + 1):
        lt = (math.lgamma(M + 1) - math.lgamma(s + 1) - math.lgamma(M - s + 1)
              + math.lgamma(n + 1) - math.lgamma(s) - math.lgamma(n - s + 2)
              + t * s * math.log((s - 1) / n))
        tot += math.exp(lt)
    return tot

print('Part A: Hall failure, t uniform candidates per row in an n-set (model of Lemma 3.4)')
print('  M   c   t     n   MC-fail/trials   exact-union-bound   lemma 1.5/n')
for (M, c, t) in [(20, 2, 5), (20, 2, 2), (20, 1.1, 3), (40, 2, 5), (40, 1.1, 2), (30, 1.5, 3)]:
    n = int(math.ceil(c * M))
    trials = 4000; fail = 0
    for _ in range(trials):
        rows = [set(random.randrange(n) for _ in range(t)) for _ in range(M)]
        if not hall_ok(rows): fail += 1
    print('%4d %4.1f %3d %5d   %5d/%d        %.3e          %.3e' % (M, c, t, n, fail, trials,
          union_bound(M, n, t), 1.5 / n))

print('\n  Sum over q in Q_large(c) = {q <= 3M^A : n_q >= cM} (Theorem 3.2 needs Hall at all of them):')
print('  M   A   c  t   sum_q 1.5/n_q (lemma)   sum_q exact-union-bound')
for (M, A, c, t) in [(50, 2, 2, 5), (50, 3, 2, 5), (100, 2, 2, 5), (100, 3, 2, 3), (100, 3, 2, 2)]:
    X = 3 * M**A
    s1 = s2 = 0.0
    # primes up to X by sieve
    sieve = bytearray([1]) * (X + 1); sieve[0] = sieve[1] = 0
    for i in range(2, int(X**0.5) + 1):
        if sieve[i]: sieve[i*i::i] = bytearray(len(sieve[i*i::i]))
    for q in range(3, X + 1, 2):
        if not sieve[q]: continue
        nq = (q - (1 if q % 4 == 1 else -1)) // 2
        if nq < c * M: continue
        s1 += 1.5 / nq
        s2 += union_bound(M, nq, t)
    print('%4d %3d %3d %2d        %.3f                 %.3e' % (M, A, c, t, s1, s2))

# ---------------- Part B -------------------------------------------------------------------------
print('\nPart B: genuine labels, Hall condition (b) of IFL at every q in Q_large(c), evidence')

def gauss_prime(p):
    for a in range(1, int(p**0.5) + 1):
        b2 = p - a * a
        b = int(round(b2**0.5))
        if b * b == b2: return (a, b)
    raise ValueError

def primroot(q):
    fac = [d for d in range(2, q) if (q - 1) % d == 0 and is_prime(d)]
    for g in range(2, q):
        if all(pow(g, (q - 1) // d, q) != 1 for d in fac): return g

def gmul(u, v, q): return ((u[0]*v[0] - u[1]*v[1]) % q, (u[0]*v[1] + u[1]*v[0]) % q)
def gpow(u, e, q):
    r = (1, 0)
    while e:
        if e & 1: r = gmul(r, u, q)
        u = gmul(u, u, q); e >>= 1
    return r

def label_table(q):
    """returns function (a,b) -> log of (a+bi)/(a-bi) in T_q (cyclic of order m), and m."""
    if q % 4 == 1:
        m = q - 1
        c = next(x for x in range(2, q) if (x * x + 1) % q == 0)
        g = primroot(q); lg = {}; y = 1
        for k in range(m): lg[y] = k; y = y * g % q
        def lab(a, b):
            x = (a + b * c) * pow((a - b * c) % q, -1, q) % q
            return lg[x]
        return lab, m
    else:
        m = q + 1
        # generator of norm-one subgroup: gamma^(q-1) for gamma generating F_{q^2}^*
        order = q * q - 1
        fac = [d for d in range(2, order + 1) if order % d == 0 and is_prime(d)]
        gam = None
        for a0 in range(q):
            for b0 in range(1, q):
                u = (a0, b0)
                if all(gpow(u, order // d, q) != (1, 0) for d in fac):
                    gam = u; break
            if gam: break
        tg = gpow(gam, q - 1, q)
        lg = {}; y = (1, 0)
        for k in range(m): lg[y] = k; y = gmul(y, tg, q)
        def lab(a, b):
            num = (a % q, b % q); den = (a % q, (-b) % q)
            nd = (den[0]*den[0] + den[1]*den[1]) % q; t = pow(nd, -1, q)
            inv = (den[0]*t % q, (-den[1])*t % q)
            return lg[gmul(num, inv, q)]
        return lab, m

M, nJ, t, c = 12, 12, 5, 2
r = 2 * M + nJ
Jc = list(range(nJ)); fcol = [nJ + 2*x for x in range(M)]; scol = [nJ + 2*x + 1 for x in range(M)]
eps = [random.choice([1, -1]) for _ in range(M)]
Aprof = [[random.randint(0, 1) if j < nJ else 0 for j in range(r)] for _ in range(M)]
for z in range(M):
    col = [random.randint(0, 1) for _ in range(M)]
    col[z] = 0 if eps[z] == 1 else 1
    for y in range(M):
        Aprof[y][scol[z]] = col[y]; Aprof[y][fcol[z]] = col[y] + (eps[z] if y == z else 0)
cand_p = [p for p in range(13, M**2 * 3) if is_prime(p) and p % 4 == 1]
primes = random.sample(cand_p, r)
pis = [gauss_prime(p) for p in primes]
n = r - M
phi = [[Aprof[y][j] for j in Jc] + [Aprof[y][scol[x]] for x in range(M)] for y in range(M)]
S = [[0] + random.sample(range(n), t - 1) for y in range(M)]   # free coordinates j_i (0 means w_0 = 0)

qmax = 1500
nfail = nq_tested = parity_bad = 0
fail_list = []
for q in range(3, qmax + 1):
    if not is_prime(q) or q in primes: continue
    m = q - (1 if q % 4 == 1 else -1)
    nq = m // 2
    lab, mm = label_table(q)
    assert mm == m
    ell = [lab(a, b) for (a, b) in pis]
    for j, p in enumerate(primes):
        nu = 1 if pow(p, (q - 1) // 2, q) == q - 1 else 0
        if ell[j] % 2 != nu: parity_bad += 1
    if nq < c * M: continue
    nq_tested += 1
    f = [ell[j] for j in Jc] + [(ell[fcol[x]] + ell[scol[x]]) % m for x in range(M)]
    kap = [sum(Aprof[y][j] * ell[j] for j in range(r)) % m for y in range(M)]
    ok = True
    for par in (0, 1):
        rows = [y for y in range(M) if kap[y] % 2 == par]
        cands = [set((kap[y] - (2 * f[S[y][i] ] if i > 0 else 0)) % m for i in range(t)) for y in rows]
        if not hall_ok(cands): ok = False
    if not ok:
        nfail += 1; fail_list.append(q)
print('M=%d r=%d primes in [13,%d], t=%d, c=%d: moduli tested (n_q >= cM, q <= %d): %d; Hall failures: %d %s'
      % (M, r, 3 * M**2, t, c, qmax, nq_tested, nfail, fail_list))
print('forced Legendre parity  ell_j(q) = [(p_j/q) = -1] mod 2  violated in %d (j,q) instances' % parity_bad)
