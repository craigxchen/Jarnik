"""Referee-2 check T5: constants of Theorem 2, the window lemma, Breusch, Prop 1.1 arithmetic, and the
prime-level bound obtained by running outside.md Theorem A.3' on the b = 8 skew-paired profile.

All comparisons of constants use 50-digit Decimal / Fraction arithmetic; the last table is arithmetic in logs.
"""
import math
from fractions import Fraction
from decimal import Decimal, getcontext
from gz import is_prime, primes_upto
getcontext().prec = 50

def asin(x):          # Decimal arcsin by series (x small)
    x = Decimal(x); s = x; term = x; k = 0
    while True:
        k += 1
        term = term * x * x * Decimal((2 * k - 1) ** 2) / Decimal((2 * k) * (2 * k + 1))
        s += term
        if abs(term) < Decimal(10) ** -45:
            return s

# (a) arcsin t <= 1.0368 t for t <= 5^(-1/2) (ratio increasing in t)
t = Decimal(1) / Decimal(5).sqrt()
print("(a) arcsin(5^-1/2)/5^-1/2 = %s  <= 1.0368: %s" % (str(asin(t) / t)[:12], asin(t) / t <= Decimal("1.0368")))
# (b) 4*1.0368*e^(1/16) <= 4.415
v = 4 * Decimal("1.0368") * (Decimal(1) / 16).exp()
print("(b) 4*1.0368*e^(1/16) = %s <= 4.415: %s" % (str(v)[:10], v <= Decimal("4.415")))
# (c) sum_{k>=2}(2k+2) x^k <= 6.3 x^2 for 0 < x <= 1/40: closed form x^2 (6/(1-x) + 2x/(1-x)^2), increasing in x
x = Fraction(1, 40)
val = x * x * (6 / (1 - x) + 2 * x / (1 - x) ** 2)
print("(c) at x=1/40: sum = %.6f x^2  (<= 6.3 x^2: %s)" % (float(val / (x * x)), val <= Fraction(63, 10) * x * x))
# (d) final constant: (4.415/C)(4x + 6.3x^2) with x <= C/40, C <= 1  ->  <= 4.415*(4 + 6.3/40)/40
d = Fraction(4415, 1000) * (4 + Fraction(63, 10) / 40) / 40
print("(d) 4.415(4 + 6.3/40)/40 = %.6f < 0.46: %s" % (float(d), d < Fraction(46, 100)))
# (e) density/window count: 2 Delta n/pi + 2 <= n + 2 needs Delta <= pi/2 (Delta <= e^-1 since W >= 4, C <= 1)
print("(e) Delta <= e^-1 < pi/2: %s" % (math.exp(-1) < math.pi / 2))
# (f) window lemma: ceil(4 r log 2) <= 2.8 r -- exact list of failures for r <= 200
fails = [r for r in range(1, 201) if math.ceil(4 * r * math.log(2)) > 2.8 * r]
print("(f) ceil(4 r log2) > 2.8 r exactly for r in", fails, "-> holds for all r >= %d; Theorem 2 has r = 9M-8 >= 64" % (max(fails) + 1))
# (g) Breusch: for every M in [8, 20000] there is a prime q = 3 mod 4, q >= 7, with M <= q+1 <= 2M
ps3 = [p for p in primes_upto(40001) if p % 4 == 3 and p >= 7]
import bisect
bad = []
for M in range(8, 20001):
    i = bisect.bisect_left(ps3, M - 1)
    if not (i < len(ps3) and ps3[i] + 1 <= 2 * M):
        bad.append(M)
print("(g) M in [8,20000] without a prime q = 3 mod 4 (q>=7) with M <= q+1 <= 2M:", bad)
# (h) Prop 1.1 arithmetic
from math import gcd
def lcm_upto(n):
    v = 1
    for k in range(1, n + 1):
        v = v * k // gcd(v, k)
    return v
print("(h) Q_2 = lcm(1..4) = %d; 12/sqrt5 = %.3f > 1; chord |(2+i)-(1+2i)| = sqrt2 = %.4f, sqrt R = 5^(1/4) = %.4f"
      % (lcm_upto(4), 12 / 5 ** 0.5, 2 ** 0.5, 5 ** 0.25))
for M in (10, 100, 1000):
    QM = lcm_upto(2 * M)
    print("    M=%d: 4(M-1) log Q_M / M^2 = %.3f  (-> 8)" % (M, 4 * (M - 1) * math.log(QM) / M ** 2))

# (i) prime-level class: outside.md A.3' run on the b = 8 skew-paired profile.
# Uniform strengthening dist(<v,phi>,(pi/4)Z) >= arcsin(Q e^{-w(v)/2}) for ALL v != 0, Q = odd lcm(<= 2M)
# (bounds every collision modulus of every monomial for ANY prime-level residue assignment, hence also the
# induced point residues, the parity law holds automatically).  Requirements:
#   characters: (4.415/C) Q [16 M^2 e^{-tau} + sum_{n>=4} (2M)^n e^{-tau n/2}(n+2)] <= 1/2   (Theorem 1, b=8)
#   completion: 1.3201 Q (exp(2r/(sqrt P - 1)) - 1) <= 1/2
# W <= r (log P + log 2 + 1/(4r)) (window lemma slack).  Report W/M^2 (limit 9*4 = 36).
def oddlcm_log(Q0):
    s = 0.0
    for q in primes_upto(Q0):
        if q == 2:
            continue
        a = 1
        while q ** (a + 1) <= Q0:
            a += 1
        s += a * math.log(q)
    return s
print("(i) prime-level class via A.3' + Theorem 1 (C = 1):")
print("      M      log Q    tau_char   tau_compl    W/M^2    W/(M log M)")
for M in (12, 44, 104, 1004, 10004, 100004):
    r = 9 * M - 8
    lQ = oddlcm_log(2 * M)
    # characters: dominant n=2 term: 4.415*Q*16 M^2 e^{-tau} <= 1/4  and the n>=4 tail is smaller once 4M^2 Q e^{-tau} < 1
    tau_c = lQ + math.log(4.415 * 16 * M * M * 4)
    # completion: exp(2r/(sqrt P -1)) - 1 <= 1/(2*1.3201 Q)  -> sqrt P - 1 >= 2r / log(1 + 1/(2.6402 Q))
    y = 1 / (2.6402 * math.exp(lQ)) if lQ < 700 else 0.0
    if y > 1e-12:
        sP = 2 * r / math.log1p(y) + 1
        tau_p = 2 * math.log(sP)
    else:
        tau_p = 2 * (math.log(2 * r * 2.6402) + lQ)
    tau = max(tau_c, tau_p)
    W = r * (tau + math.log(2) + 1 / (4 * r))
    print("  %7d  %9.1f  %9.1f   %9.1f   %7.2f   %10.1f" % (M, lQ, tau_c, tau_p, W / M ** 2, W / (M * math.log(M))))
print("T5 DONE")
