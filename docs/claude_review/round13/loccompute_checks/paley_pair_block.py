# two-cut block cut out by two Paley columns (support 9 = (3,3,3) Venn), L = (2,2), degrees <= 8: J vs J_w(R=1)
import time, itertools
from locrun import *
t0 = time.time()
d = 9; S1 = list(range(0, 6)); S2 = list(range(3, 9))
A = models.box(d, *models.cuts(d, [S1, S2], [2, 2])[1:])
m = 8; ring = Ring(d, m)
H = hilbert_modp(A, ring, m); JA = [ring.N[i]-H[i] for i in range(m+1)]
gw = width_gens(ring, width_dirs(A, d, 1, m))
Jw = ideal_dims_modp(ring, gw, m)
print(f"two-cut (3,3,3) L=(2,2) d=9 |A|={len(A)} H={H}\n  J(A)={JA}\n  J_w ={Jw} cum_eq={cum_eq(Jw, JA)} [{time.time()-t0:.0f}s]")
