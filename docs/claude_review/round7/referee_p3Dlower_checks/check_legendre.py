"""Referee check: parity of log(z/conj z) in T_q equals [(Norm z / q) = -1] (sharp.md Prop 2.1(3)), and
i is a square in T_q iff q = +-1 mod 8, for all odd primes q < 200 and z = a+bi, 0<=a<25, 1<=b<25, q not | Norm z."""
from sympy import primerange, legendre_symbol
from tq import T
ok = True; n = 0
for q in primerange(3, 200):
    t = T(q)
    ok &= ((t.iota % 2 == 0) == (q % 8 in (1, 7)))
    for a in range(25):
        for b in range(1, 25):
            N = a * a + b * b
            if N % q == 0: continue
            l = t.L(a, b); n += 1
            ok &= ((l % 2) == (1 if legendre_symbol(N % q, q) == -1 else 0))
print("parity of log(z/conj z) == [(Norm z/q) = -1], and iota even iff q = +-1 mod 8:", ok, "cases", n)
print("ALL OK" if ok else "FAIL")
