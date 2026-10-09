"""Referee check R6: Theorem 2 of construct.md -- the constant chain of steps 1-4 and an independent instance at
the smallest admissible order M = 8 (q = 7), C = 1, Q0 = 2M, P = P_* exactly as in the theorem.

Constant chain (50-digit Decimal):
  (a) arcsin t <= 1.0368 t for 0 <= t <= 5^(-1/2) (max of arcsin(t)/t is at the endpoint, arcsin convex);
  (b) 4 * 1.0368 * e^(1/16) <= 4.415;
  (c) sum_{n>=4 even} (n+2) x^(n/2) <= 6.3 x^2 for 0 < x <= 1/40;
  (d) (4.415)(4 + 6.3/40)/40 < 0.46.
Instance (exact primes, Decimal logs): L = 8th prime = 1 mod 4 above 16, P = ceil(max(160 M^2 L / C, 21 r^2)),
r = 64 primes = 1 mod 4 in [P, P e^(1/(4r))]; then the union bound of step 4 with the proved Theorem 1 margins
(g(2) = 4, g(n) = 2n for n >= 4) and the count (2M)^n, the completion constant 1.3201 (Pi - 1) < 1, and
W/(M log M)."""
from decimal import Decimal, getcontext
import math
getcontext().prec = 50

def dexp(x): return x.exp()
# (a) arcsin(t)/t increasing on (0,1): check value at 5^-1/2 by series
t = Decimal(1) / Decimal(5).sqrt()
s = Decimal(0); term = t; k = 0
while True:     # arcsin t = sum (2k)!/(4^k (k!)^2 (2k+1)) t^(2k+1)
    c = Decimal(math.comb(2*k, k)) / Decimal(4**k) / Decimal(2*k + 1)
    tt = c * t**(2*k + 1)
    s += tt
    if tt < Decimal(10)**-45: break
    k += 1
print("(a) arcsin(5^-1/2)/5^-1/2 = %s  (<= 1.0368: %s)" % (str(s / t)[:12], s / t <= Decimal("1.0368")))
assert s / t <= Decimal("1.0368")
v = 4 * Decimal("1.0368") * dexp(Decimal(1) / 16)
print("(b) 4*1.0368*e^(1/16) = %s  (<= 4.415: %s)" % (str(v)[:10], v <= Decimal("4.415")))
assert v <= Decimal("4.415")
x = Decimal(1) / 40
ser = sum((2*k + 2) * x**k for k in range(2, 200))
print("(c) sum_{k>=2} (2k+2) x^k at x=1/40: %s = %s x^2 (<= 6.3 x^2: %s); ratio is increasing in x"
      % (str(ser)[:12], str(ser / x**2)[:8], ser <= Decimal("6.3") * x**2))
assert ser <= Decimal("6.3") * x**2
d = Decimal("4.415") * (4 + Decimal("6.3") / 40) / 40
print("(d) 4.415*(4+6.3/40)/40 = %s < 0.46: %s" % (str(d)[:8], d < Decimal("0.46")))
assert d < Decimal("0.46")

# ---- instance M = 8
def is_prime(n):
    if n < 2: return False
    for p in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
        if n % p == 0: return n == p
    d_, s_ = n - 1, 0
    while d_ % 2 == 0: d_ //= 2; s_ += 1
    for a in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
        y = pow(a, d_, n)
        if y in (1, n - 1): continue
        for _ in range(s_ - 1):
            y = y * y % n
            if y == n - 1: break
        else: return False
    return True

M = 8; b = 8; r = b*(M - 1) + M; C = Decimal(1); Q0 = 2*M
ells = []; t_ = Q0 + 1
while len(ells) < M:
    if t_ % 4 == 1 and is_prime(t_): ells.append(t_)
    t_ += 1
L = ells[-1]
Pstar = max(-(-160 * M * M * L // 1), 21 * r * r)
P = Pstar
# slide the window start until r primes = 1 mod 4 fit in [P, P e^(1/(4r))]
while True:
    hi = Decimal(P) * dexp(Decimal(1) / (4 * r))
    ps = [n for n in range(P, int(hi) + 1) if n % 4 == 1 and is_prime(n)]
    if len(ps) >= r: ps = ps[:r]; break
    P += 1
w = [Decimal(p).ln() for p in ps]
tau = Decimal(P).ln(); eta = w[-1] - tau
assert eta <= Decimal(1) / (4 * r)
W = sum(w)
Pi = Decimal(1)
for p in ps: Pi *= 1 + 2 / (Decimal(p).sqrt() - 1)
print("instance M=%d r=%d L=%d P_*=%d window start P=%d, primes %d..%d, eta=%.2e <= 1/(4r)=%.2e"
      % (M, r, L, Pstar, P, ps[0], ps[-1], eta, 1 / (4 * r)))
print("  completion: 1.3201 (Pi - 1) = %s < 1" % str(Decimal("1.3201") * (Pi - 1))[:8])
assert Decimal("1.3201") * (Pi - 1) < 1
xx = 4 * M * M * Decimal(L) / Decimal(P)
assert xx <= C / 40
# union bound: per character (4*1.0368/C) e^{eta r/4} L^{n/2} e^{-tau g/4} (n+2); count (2M)^n
tot = Decimal(0)
for n in range(2, 2000, 2):
    g = 4 if n == 2 else 2 * n
    term = (4 * Decimal("1.0368") / C) * dexp(eta * r / 4) * Decimal(2*M)**n * Decimal(L)**n / Decimal(L)**(n // 2) \
        if False else (4 * Decimal("1.0368") / C) * dexp(eta * r / 4 + n * (Decimal(2*M).ln() + Decimal(L).ln() / 2) - tau * g / 4) * (n + 2)
    tot += term
print("  x = 4M^2 L/P = %s (<= 1/40); union bound P(fail) <= %s (analytic bound 0.46)" % (str(xx)[:8], str(tot)[:8]))
assert tot < Decimal("0.46")
print("  W = %.2f, W/(M log M) = %.2f;  3 M log M = %.1f <= W + 14 M (two_thirds): %s"
      % (W, float(W) / (M * math.log(M)), 3 * M * math.log(M), 3 * M * math.log(M) <= float(W) + 14 * M))
# asymptotic constant: log P_* / log M with L ~ M-th prime = 1 mod 4 above 2M
for Mq in (12, 104, 1004, 10004):
    ls = []; t_ = 2*Mq + 1; cnt = 0
    while cnt < Mq:
        if t_ % 4 == 1 and is_prime(t_): cnt += 1; last = t_
        t_ += 1
    rr = 9*Mq - 8
    lp = math.log(max(160 * Mq * Mq * last, 21 * rr * rr))
    print("  M=%6d: L=%d, log P_*/log M = %.3f, r log P_*/(M log M) = %.2f (limit 27)" % (Mq, last, lp / math.log(Mq), rr * lp / (Mq * math.log(Mq))))
print("R6 DONE")
