"""W/(M log M) of Theorem 5.1's fake circles from the explicit P_* (C = 1), for large M:
   W <= r (log P_* + 1/(4r)),  r = 2(M-1),  P_* = max(160000 M^4 L^2, 21 r^2),
   L = M-th prime = 1 mod 4 above Q0 = 2M (A_1), or L = 1 (A_0).  L computed exactly by a sieve.
Floating point (double) in the logarithms; this is arithmetic on the proved bound, not a proof step."""
import math
import numpy as np


def primes_1mod4_above(Q0, count):
    N = max(100, int((Q0 + 4 * count * math.log(max(count, 3)) * 2) * 1.5))
    while True:
        s = np.ones(N + 1, dtype=bool); s[:2] = False
        for i in range(2, int(N ** 0.5) + 1):
            if s[i]:
                s[i * i::i] = False
        ps = [p for p in np.nonzero(s)[0] if p > Q0 and p % 4 == 1]
        if len(ps) >= count:
            return int(ps[count - 1])
        N *= 2


print(' M        A_0: W/(MlogM)   A_1(2M): W/(MlogM)   L')
for M in [8, 12, 20, 44, 104, 1004, 10004, 100004, 1000004]:
    r = 2 * (M - 1)
    L = primes_1mod4_above(2 * M, M)
    for lab, LL in (('A0', 1), ('A1', L)):
        P = max(160000 * M ** 4 * LL ** 2, 21 * r * r)
        W = r * (math.log(P) + 1 / (4 * r))
        if lab == 'A0':
            a0 = W / (M * math.log(M))
        else:
            a1 = W / (M * math.log(M))
    print(f'{M:8d}   {a0:8.2f}            {a1:8.2f}          {L}')
