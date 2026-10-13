"""Referee checks (independent):
G. Prop 3.1 scope: additive labels ell = g_q*lambda + iota*u are parity-respecting at q = +-1 mod 8 only if
   nu(q) in {0, lambda mod 2}. For random genuine primes, count the q = +-1 mod 8 <= X where this can hold for the best lambda.
H. Thm 3.2 evidence: least prime q = +-1 mod 8, q not | n, with (n/q) = -1 for n = p*Norm(zeta) not a square / twice a square.
I. Cor 4.4 algebra (sympy)."""
import random, math, itertools
from sympy import primerange, legendre_symbol, isprime, symbols, solve, simplify, Rational
random.seed(3)
# G
for (r, X, lo) in [(20, 300, 10**6), (60, 1000, 10**7)]:
    P = set()
    while len(P) < r:
        c = random.randrange(lo, 2 * lo)
        if c % 4 == 1 and isprime(c): P.add(c)
    P = sorted(P)
    Q1 = [q for q in primerange(3, X + 1) if q % 8 in (1, 7)]
    Q3 = [q for q in primerange(3, X + 1) if q % 8 in (3, 5)]
    nus = {q: tuple(1 if legendre_symbol(p % q, q) == -1 else 0 for p in P) for q in Q1}
    zero = tuple([0] * r)
    best = max(sum(1 for q in Q1 if nus[q] in (zero, lb)) for lb in set(nus.values()) | {zero})
    th1 = sum(math.log(q) for q in Q1); th3 = sum(math.log(q) for q in Q3)
    print(f"G: r={r}, X={X}: #q=+-1 mod 8: {len(Q1)}, distinct nu(q): {len(set(nus.values()))}; "
          f"max #q=+-1 mod 8 compatible with ONE additive lambda: {best}; log-prod over q=+-3 mod 8: {th3:.1f} vs all odd: {th1+th3:.1f}")
    # pairs (j,k) with (p_j p_k / q) = 1 for all q = +-1 mod 8 <= X
    bad = sum(1 for j, k in itertools.combinations(range(r), 2) if all(nus[q][j] == nus[q][k] for q in Q1))
    print(f"   pairs j<k with equal Legendre symbols at all q=+-1 mod 8 <= X: {bad}")
# H
worst = 0; wcase = None; cnt = 0
for p in [p for p in primerange(5, 3000) if p % 4 == 1]:
    for a in range(0, 60):
        for b in range(1, 60):
            nz = a * a + b * b
            if nz >= p: continue
            n = p * nz
            s = math.isqrt(n); s2 = math.isqrt(n // 2)
            if s * s == n or (n % 2 == 0 and 2 * s2 * s2 == n): continue
            q = next(q for q in primerange(3, 10**6) if q % 8 in (1, 7) and n % q and legendre_symbol(n % q, q) == -1)
            cnt += 1
            ratio = q / math.log(n) ** 2
            if ratio > worst: worst, wcase = ratio, (p, (a, b), n, q)
print(f"H: {cnt} pairs (p, zeta) with Norm zeta < p: max over them of (least q=+-1 mod 8 with (n/q)=-1)/(log n)^2 = {worst:.3f} at {wcase}")
# I
k, T, h, lM, lam, tau = symbols('k T h lM lam tau', positive=True)
cond = (2*T*k + T + 4)*(h*lM + Rational(7, 10)) + T*(k + 1)*lam - tau
k1 = solve(cond, k)[0]
target = (tau - (T + 4)*(h*lM + Rational(7, 10)) - T*lam) / (T*(2*h*lM + Rational(14, 10) + lam))
print("I: Cor 4.4 k_1 matches (4.3):", simplify(k1 - target) == 0)
