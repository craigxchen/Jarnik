"""Tight weight-<=2 part of the completion sum: sum over tight (v,s) of the excluded slab measure 2d/pi."""
import sys, math, random
from rcore import rand_primes_like_author, gauss, Lprime, tight_list
for (M, b, K) in [(164, 8, 4), (644, 1, 4)]:
    P = int(M ** K); ps = rand_primes_like_author(b * (M - 1), P, 2 * P, random.Random(0))
    gp = [gauss(p) for p in ps]; L = Lprime(2 * M, set(ps))
    rels = tight_list(gp, L)
    tot = sum(2 * math.asin(math.sqrt(float(s2))) / math.pi for v, s, sg, s2 in rels)
    pairs = sum(2 * math.asin(math.sqrt(float(s2))) / math.pi for v, s, sg, s2 in rels if len(v) == 2)
    print(f'M={M} b={b} K={K}: tight={len(rels)}  sum 2d/pi = {tot:.1f}  (pairs only {pairs:.1f})')
