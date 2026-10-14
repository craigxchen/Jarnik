# Sylvester profile n=2^k at m'=ceil(4w): is the 4-column span cheap?  critical transversal order for Prop 5.1.
import math
from fractions import Fraction as F
from params import j_mprime
for k in range(12, 21):
    n = 2**k; d = n-1; r = n-1
    w, j, mp_ = j_mprime(n, r)
    w = float(w)
    out = []
    for c in (3, 4):
        cost = 2*j*c
        cheap = cost <= mp_-1
        P = 2**c * (n//2**c)*(n//2**c - 1)//2
        N2 = n*(n-1)//2
        ocrit = F(d*P, N2)   # Prop 5.1 needs o <= d P/N2
        out.append(f"{c}col: cost={cost} cheap={cheap} o_lower={mp_-cost} o_crit={float(ocrit):.1f}")
    print(f"k={k} n={n} j={j} m'={mp_} | " + " | ".join(out))
