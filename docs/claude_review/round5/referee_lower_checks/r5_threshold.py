"""Referee r5 (illustration, floating point): the explicit union bound written in the proof of
lower.md Theorem 5.1(a) (C = 1), as a function of L = log M, at the worst admissible tau = log Z,
Z = M^6 (log M)^6, Y = 3 M log M:
  pairs   (12 pi/C) e^(1/4) * 81 / L^2
  ternary (8 pi/C) e^(1/4) * sum_{h>=4} (1+h) x^h,  x = 18/L
(the amplitude term is negligible).  Prints the smallest L (step 1) where the written total is < 1,
and W/(M log M) = 8 * (6 + 6 log L / L) at that point (the o(1) of '48 + o(1)')."""
import math
def total(L, C=1.0):
    e = math.exp(0.25)
    pairs = 12 * math.pi / C * e * 81 / L ** 2
    x = 18 / L
    tern = 8 * math.pi / C * e * sum((1 + h) * x ** h for h in range(4, 400)) if x < 1 else float('inf')
    return pairs, tern
for L in range(20, 200):
    p, t = total(L)
    if p + t < 1:
        print("C=1: written bound < 1 first at log M = %d (M ~ 10^%.0f): pairs %.3f, ternary %.3f;"
              " W/(M log M) there ~ %.1f" % (L, L / math.log(10), p, t, 8 * (6 + 6 * math.log(L) / L)))
        break
for C in (0.1, 0.01):
    for L in range(20, 2000):
        p, t = total(L, C)
        if p + t < 1:
            print("C=%g: first log M = %d" % (C, L)); break
