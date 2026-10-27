"""Sanity for slab_iv: with the genuine label at q0, a tiny box around phi = arg pi_j is never discarded."""
import sys, math
from fractions import Fraction as Fr
from slab_iv import constraints, check_box
from sllib import two_sq, genuine
Fpath, q0, K = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]); ps = [int(x) for x in sys.argv[4:]]
m, L = genuine(ps, q0)
cons = constraints(Fpath, q0, K, L, ps)
t0 = [Fr(math.atan2(b, a)/math.pi).limit_denominator(10**15) for a, b in map(two_sq, ps)]
eps = Fr(1, 10**12)
print("genuine label", L, "tiny box at genuine angles discarded by constraint:", check_box(cons, [t-eps for t in t0], [t+eps for t in t0]))
