"""Exact pair slacks sigma_xy = sum_S mu_S (|lambda(S)| - 1/2) for lambda = e_x - e_y, sdc_q (mu = 1),
and r = number of patterns; also 4- and 6-term balanced lambda statistics (sample) for the light/heavy split."""
import itertools, random
from fractions import Fraction
from flib import sdc_patterns
for q in (7, 11, 19, 23, 31, 43):
    M, pats = sdc_patterns(q)
    r = len(pats)
    Ss = [set(P) for P in pats]
    sig = []
    for x, y in itertools.combinations(range(M), 2):
        s = sum(Fraction(abs((x in S) - (y in S))) - Fraction(1, 2) for S in Ss)
        sig.append(s)
    random.seed(1)
    def lamsig(pos, neg):
        return sum(Fraction(abs(sum(p in S for p in pos) - sum(t in S for t in neg))) - Fraction(1, 2) for S in Ss)
    s4 = [lamsig(*[random.sample(range(M), 4)[i::2] for i in range(2)]) for _ in range(200)]
    s6 = [lamsig(*(lambda v: (v[:3], v[3:]))(random.sample(range(M), 6))) for _ in range(200)]
    s8 = [lamsig(*(lambda v: (v[:4], v[4:]))(random.sample(range(M), 8))) for _ in range(200)] if M >= 8 else [0]
    print(f"sdc{q}: n={M} r={r} pair sigma: min {min(sig)} max {max(sig)} | 4-term sigma/r mean {float(sum(s4))/len(s4)/r:.3f} "
          f"| 6-term {float(sum(s6))/len(s6)/r:.3f} | 8-term {float(sum(s8))/len(s8)/r:.3f}  (light iff < 0.5)")
