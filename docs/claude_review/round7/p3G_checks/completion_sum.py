"""Exact lower bound for the weight-<=2 part of the completion sum on two verified prime sets:
sum over tight relations (v,s) of (2/pi) arcsin(nu_s m_(v,s) e^(-w/2)) = (2/pi) d_e  (Lemma 2.3)."""
import random, math
from fractions import Fraction
from gtight import two_squares, tight_relations, S_value
from trend import rand_primes
for (M, b, K, seed) in [(164, 8, 4, 0), (644, 1, 4, 0)]:
    rng = random.Random(seed)
    r = b * (M - 1); P = int(M ** K)
    ps = rand_primes(r, P, 2 * P, rng)
    gp = []
    for p in ps:
        a, bb = two_squares(p)
        if rng.random() < 0.5: bb = -bb
        gp.append((a, bb))
    rels = tight_relations(gp, 2 * M)
    tot = sum((2 / math.pi) * math.asin(math.sqrt(float(S_value(x)))) for x in rels)
    totP = sum((2 / math.pi) * math.asin(math.sqrt(float(S_value(x)))) for x in rels if x[0] == 'P')
    print(f'M={M} b={b} K={K}: tight relations {len(rels)}, sum of bad-slab measures over tight weight<=2 relations = {tot:.1f} (pairs alone {totP:.1f})')
