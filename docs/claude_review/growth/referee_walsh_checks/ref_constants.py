"""Referee check 4: integer constants of Cor B.2 and Thm B''."""
import math
from math import comb, isqrt, log
def s(M):  # floor((1+sqrt(4M-3))/2)+8, exact
    r = isqrt(4*M-3)
    v = (1 + r)//2
    # floor((1+sqrt(x))/2): adjust in case sqrt non-integer
    while (2*(v+1)-1)**2 <= 4*M-3: v += 1
    while (2*v-1)**2 > 4*M-3: v -= 1
    return v + 8
orders = [M for M in range(4, 5000, 4)]
first_fail = [M for M in orders if not (M/2 > s(M))]
print("Hadamard orders (multiples of 4) with M/2 <= s_M:", first_fail)
# 14 smallest split primes
sp = [p for p in range(5, 200) if p % 4 == 1 and all(p % d for d in range(2, int(p**0.5)+1))][:14]
print("14 smallest split primes:", sp, " sum log =", round(sum(map(log, sp)), 3), " 36 log 2 =", round(36*log(2), 3))
# Theorem B'' pigeonhole for all multiples of 4 from 4096 to 70000 (step 4), with conservative counts
bad = []
for M in range(4096, 70001, 4):
    f = int(M/(40*log(M)))
    L = M//60
    heavy = (32*L)//7      # heavy < 32L/7  => heavy <= floor(32L/7) (if 32L/7 integer, strict gives -1)
    A = -(-(M-f)//3)       # |A| >= (M-f)/3
    n = A - heavy
    # log binom(n,L) via lgamma vs log((2L+1)^f M)
    lb = math.lgamma(n+1) - math.lgamma(L+1) - math.lgamma(n-L+1)
    rb = f*log(2*L+1) + log(M)
    if not lb > rb + 1: bad.append(M)
print("Thm B'' pigeonhole failures in [4096,70000]:", bad[:5], len(bad))
# smallest M where it would fail (to see the margin)
for M in [256, 512, 1024, 2048, 3000, 4096]:
    f = int(M/(40*log(M))); L = M//60; heavy=(32*L)//7; A=-(-(M-f)//3); n=A-heavy
    lb = math.lgamma(n+1)-math.lgamma(L+1)-math.lgamma(n-L+1); rb=f*log(2*L+1)+log(M)
    print("  M=%d f=%d L=%d n=%d  log binom=%.1f  log records=%.1f" % (M,f,L,n,lb,rb))
