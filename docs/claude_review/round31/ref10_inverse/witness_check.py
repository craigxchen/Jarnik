#!/usr/bin/env python3
"""Exact Gaussian-integer check of X.3 on actual IX.3-type witnesses built from a minimal V.
(1) global pair identities A/conj(A) in units, for all pairs of owners;
(2) log H_i / T for owners, and the bound of case 2;
(3) case 1: twist one witness by an exterior split prime beta (norm l0 > X0); find the least
    exterior odd prime l at which some pair congruence fails; check sum_{l < l_fail} log l <= log(2|A|).
"""
import sys, math, random
sys.path.insert(0, '.')
from minV import setup, solve_col
from sympy import isprime, primerange

def gmul(a, b): return (a[0]*b[0]-a[1]*b[1], a[0]*b[1]+a[1]*b[0])
def gconj(a): return (a[0], -a[1])
def gpow(a, e):
    r = (1, 0)
    while e:
        if e & 1: r = gmul(r, a)
        a = gmul(a, a); e >>= 1
    return r
def gnorm(a): return a[0]*a[0]+a[1]*a[1]
UNITS = [(1, 0), (0, 1), (-1, 0), (0, -1)]
def gprime(p):
    for x in range(1, int(math.isqrt(p)) + 1):
        y2 = p - x*x; y = math.isqrt(y2)
        if y*y == y2: return (x, y)
    raise ValueError

def main(q):
    I, K, S, a, D, R = setup(q)
    n = len(I)
    V = [[0]*n for _ in range(n)]
    for c in range(n):
        rhs = [R[r][c] for r in range(len(D))]
        val, x = solve_col(D, rhs, n, list(range(n)))
        for r in range(n): V[r][c] = x[r]
    P = q**3
    ps = [p for p in primerange(P, 2*P) if p % 4 == 1][:n]
    pi = [gprime(p) for p in ps]
    def G(k, extra=None):
        g = (1, 0)
        for l in range(n):
            v = V[k][l]
            if v > 0: g = gmul(g, gpow(pi[l], v))
            elif v < 0: g = gmul(g, gpow(gconj(pi[l]), -v))
        if extra is not None: g = gmul(g, extra)
        return g
    Gs = [G(k) for k in range(n)]
    T = sum(math.log(p) for p in ps)
    def Apair(i, ip, Glist):
        ii, iip = I.index(i), I.index(ip)
        A = (1, 0)
        for k in range(n):
            dk = K[(i, I[k])] - K[(ip, I[k])]
            if dk > 0: A = gmul(A, gpow(gmul(pi[k], gmul(Glist[k], Glist[k])), dk))
            elif dk < 0: A = gmul(A, gpow(gmul(gconj(pi[k]), gmul(gconj(Glist[k]), gconj(Glist[k]))), -dk))
        return gmul(A, gmul(gconj(pi[ii]), pi[iip]))
    ok = 0; npair = 0
    for x in range(len(S)):
        for y in range(x+1, len(S)):
            A = Apair(S[x], S[y], Gs); Ab = gconj(A); npair += 1
            if any(gmul(u, Ab) == A for u in UNITS): ok += 1
    logH = sorted(math.log(gnorm(Gs[I.index(i)])) / T for i in S)
    print(dict(q=q, pairs=npair, exact_identities=ok, owner_logH_over_T=[round(v, 3) for v in logH]))
    # case 1
    random.seed(q)
    l0 = next(p for p in primerange(10**6, 10**7) if p % 4 == 1 and p not in ps)
    beta = gprime(l0)
    kstar = I.index(S[0])
    G2 = list(Gs); G2[kstar] = G(kstar, beta)
    worst = None
    for x in range(len(S)):
        for y in range(x+1, len(S)):
            i, ip = S[x], S[y]
            A = Apair(i, ip, G2); Ab = gconj(A)
            if any(gmul(u, Ab) == A for u in UNITS):
                continue
            # for each FIXED unit u: least exterior odd prime l with A != u conj(A) mod l
            L = -1.0; lfail = None
            for u in UNITS:
                w = gmul(u, Ab); Lu = 0.0
                for l in primerange(3, 10**6):
                    if l in ps or l == l0: continue
                    if not ((A[0]-w[0]) % l == 0 and (A[1]-w[1]) % l == 0):
                        break
                    Lu += math.log(l)
                if Lu > L: L, lfail = Lu, l
            bound = math.log(2) + 0.5 * math.log(gnorm(A))
            worst = (i, ip, lfail, round(L, 2), round(bound, 2)) if worst is None or L > worst[3] else worst
            assert L <= bound
    print(dict(q=q, case1_twisted_label=S[0], worst_pair=worst))

if __name__ == '__main__':
    for q in map(int, sys.argv[1:]): main(q)
