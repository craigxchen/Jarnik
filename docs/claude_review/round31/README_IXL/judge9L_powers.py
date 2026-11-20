# Judge check (ii): bounded-exponent Gaussian extraction (Part IX.0, "bounded powers on Gaussian circles").
# Own code. Exact Gaussian-integer arithmetic (pairs of Python ints), exact Fractions for the probability inequalities.
import math, random, itertools
from fractions import Fraction
from math import log
random.seed(9112026)
out = []
def say(*a):
    s = " ".join(str(x) for x in a); print(s); out.append(s)
def gmul(x, y): return (x[0]*y[0]-x[1]*y[1], x[0]*y[1]+x[1]*y[0])
def gconj(x): return (x[0], -x[1])
def gpow(x, e):
    r = (1, 0)
    for _ in range(e): r = gmul(r, x)
    return r
def gdiv(x, y):
    n = y[0]**2+y[1]**2; m = gmul(x, gconj(y))
    if m[0] % n or m[1] % n: return None
    return (m[0]//n, m[1]//n)
def gsub(x, y): return (x[0]-y[0], x[1]-y[1])
def gadd(x, y): return (x[0]+y[0], x[1]+y[1])
def sc(c, x): return (c*x[0], c*x[1])
def primes(n):
    return [p for p in range(2, n) if all(p % d for d in range(2, int(p**.5)+1))]
SPLIT = [p for p in primes(200) if p % 4 == 1]
def gen(p):
    for a in range(1, p):
        b = math.isqrt(p-a*a)
        if b*b == p-a*a and a > b: return (a, b)
def factor_support(n, ps):
    return {p for p in ps if n % p == 0}
def vpi(x, pi, p):
    # valuation of Gaussian x at the prime pi (above p)
    v = 0
    while True:
        y = gdiv(x, pi)
        if y is None: return v
        x = y; v += 1

stats = dict(configs=0, pairs=0, pair_fail=0, ineq_fail=0, triples=0, hidden_triples=0, case_fail=0, recover_fail=0)
for trial in range(120):
    E = random.choice([2, 3])
    s = random.randint(2, 4)
    ps = random.sample(SPLIT, s); pis = {p: gen(p) for p in ps}
    e = {p: random.randint(1, E) for p in ps}
    W = sum(e[p]*log(p) for p in ps)
    eta = random.choice([(1,0),(0,1),(-1,0),(0,-1)])
    allocs = list(itertools.product(*[range(e[p]+1) for p in ps]))
    b = min(len(allocs), random.randint(4, 9))
    rows = random.sample(allocs, b)
    stats['configs'] += 1
    def z(a):
        r = eta
        for p, ap in zip(ps, a): r = gmul(r, gmul(gpow(pis[p], ap), gpow(gconj(pis[p]), e[p]-ap)))
        return r
    Z = [z(a) for a in rows]
    def H(ai, aj):
        r = (1, 0)
        for p, x, y in zip(ps, ai, aj):
            r = gmul(r, gmul(gpow(pis[p], min(x, y)), gpow(gconj(pis[p]), e[p]-max(x, y))))
        return r
    # pair factorisation z_j - z_i = 2 i eta H t, t in Z\{0}; |H|^2 = exp(W - D)
    T = {}
    for i in range(b):
        for j in range(b):
            if i == j: continue
            stats['pairs'] += 1
            h = H(rows[i], rows[j])
            q = gdiv(gsub(Z[j], Z[i]), gmul((0, 2), gmul(eta, h)))
            D = sum(abs(x-y)*log(p) for p, x, y in zip(ps, rows[i], rows[j]))
            if q is None or q[1] != 0 or q[0] == 0 or abs(log(h[0]**2+h[1]**2) - (W-D)) > 1e-9:
                stats['pair_fail'] += 1
            else:
                T[(i, j)] = q[0]
    # exact inequalities at each prime: P(a1=a2!=a3) >= delta^2/e_p, and E|a1-a2| <= e_p delta
    for k, p in enumerate(ps):
        cnt = [sum(1 for a in rows if a[k] == v) for v in range(e[p]+1)]
        al = [Fraction(c, b) for c in cnt]
        delta = sum(x*(1-x) for x in al)
        P = sum(x*x*(1-x) for x in al)
        Eabs = sum(al[u]*al[v]*abs(u-v) for u in range(e[p]+1) for v in range(e[p]+1))
        # direct count of P(a1=a2!=a3) with replacement
        Pd = Fraction(sum(1 for i in range(b) for j in range(b) for l in range(b)
                          if rows[i][k] == rows[j][k] != rows[l][k]), b**3)
        if not (P == Pd and P*e[p] >= delta*delta and Eabs <= e[p]*delta and P*E >= delta*delta):
            stats['ineq_fail'] += 1
    # triples: division by componentwise gcd, hidden slots, recovery of the hidden set from the equation
    for (i1, i2, i3) in itertools.permutations(range(b), 3):
        stats['triples'] += 1
        a1, a2, a3 = rows[i1], rows[i2], rows[i3]
        G3 = (1, 0)
        for p, x, y, w in zip(ps, a1, a2, a3):
            G3 = gmul(G3, gmul(gpow(pis[p], min(x, y, w)), gpow(gconj(pis[p]), e[p]-max(x, y, w))))
        R12 = gdiv(H(a1, a2), G3); R23 = gdiv(H(a2, a3), G3); R13 = gdiv(H(a1, a3), G3)
        hidden = [p for p, x, y, w in zip(ps, a1, a2, a3) if x == y != w]
        # case check: at hidden p, R12 carries exactly conj(pi)^{w-x} (w>x) or pi^{x-w} (w<x); R13,R23 nothing at p
        for p, x, y, w in zip(ps, a1, a2, a3):
            if x == y != w:
                pi, pb = pis[p], gconj(pis[p])
                want = (0, w-x) if w > x else (x-w, 0)
                got = (vpi(R12, pi, p), vpi(R12, pb, p))
                if got != want or vpi(R13, pi, p) or vpi(R13, pb, p) or vpi(R23, pi, p) or vpi(R23, pb, p):
                    stats['case_fail'] += 1
        # triangle identity t12 R12 + t23 R23 = t13 R13
        t12, t23, t13 = T[(i1, i2)], T[(i2, i3)], T[(i1, i3)]
        if gadd(sc(t12, R12), sc(t23, R23)) != sc(t13, R13): stats['recover_fail'] += 1; continue
        if not hidden: continue
        stats['hidden_triples'] += 1
        # exposed part U of R12: product over non-hidden primes of its factors (known from labels)
        U = (1, 0)
        for p, x, y, w in zip(ps, a1, a2, a3):
            if p in hidden: continue
            m, M = min(x, y, w), max(x, y, w)
            U = gmul(U, gmul(gpow(pis[p], min(x, y)-m), gpow(gconj(pis[p]), M-max(x, y))))
        X = gdiv(gsub(sc(t13, R13), sc(t23, R23)), sc(t12, U))
        if X is None: stats['recover_fail'] += 1; continue
        n = X[0]**2 + X[1]**2
        if factor_support(n, ps) != set(hidden) or math.prod(p**abs(x-w) for p, x, y, w in zip(ps, a1, a2, a3) if p in hidden) != n:
            stats['recover_fail'] += 1
say("(ii) bounded-exponent Gaussian extraction:", stats)
open(__file__.replace('.py', '.out'), 'w').write("\n".join(out) + "\n")
