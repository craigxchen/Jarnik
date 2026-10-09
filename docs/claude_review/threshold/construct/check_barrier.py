"""Independent spot-check of the slice bound reg Q^M(p,n) <= phi_M(p,n) (upper.md Thm 3.1)
and of max reg/t = 2 - 1/ceil(M/2).  Mod-q rank gives rigorous UPPER bounds on reg."""
from fractions import Fraction as Fr
from clib import reg2, proj
from gens import Qslice
import math, sys
def phi(M,p,n):
    return min(Fr(M-1)*(Fr(p,M-j)+Fr(n,j)) for j in range(1,M))
tmax = {3:10, 4:8, 5:6, 6:5}
for M in (3,4,5,6):
    best = Fr(0)
    for t in range(1, tmax[M]+1):
        for n in range(0, t+1):
            p = t-n
            S = Qslice(M,p,n)
            r = reg2(proj(S))
            ph = phi(M,p,n)
            assert r <= ph, (M,p,n,r,ph)
            best = max(best, Fr(r,t))
            print(f"M={M} p={p} n={n} |Q|={len(S)} reg={r} phi={float(ph):.3f}", flush=True)
    k = math.ceil(M/2)
    print(f"M={M}: max reg/t over computed = {best} ; 2-1/ceil(M/2) = {2-Fr(1,k)}", flush=True)
