import sys, time, itertools
import numpy as np
from locrun import *
k, L, m = 3, int(sys.argv[1]), int(sys.argv[2])
d, cols, Ls = models.sylvester(k, L)
A = models.box(d, cols, Ls)
ring = Ring(d, m)
H = hilbert_modp(A, ring, m); JA = [ring.N[i]-H[i] for i in range(m+1)]
print(f"Syl(8,{L}) |A|={len(A)} m={m} H={H} J(A)={JA}")
gw = width_gens(ring, width_dirs(A, d, 2, m))
# lines {a,b,a+b}
lines = sorted(set(tuple(sorted((a,b,a^b))) for a in range(1,8) for b in range(1,8) if a!=b))
half = []
for ln in lines:
    v = (np.array(cols[ln[0]-1])+np.array(cols[ln[1]-1])+np.array(cols[ln[2]-1]))
    assert (v%2==0).all(); half.append(list(v//2))
print("half-sum vectors:", half[:3], "widths:", [len(image(A,[h])) for h in half])
cache={}
for name, fams in [("lines (3 cols)", [[cols[a-1] for a in ln] for ln in lines]),
                   ("pairs (col, half-sum of a line through it)", [[cols[a-1], half[i]] for i,ln in enumerate(lines) for a in ln]),
                   ("pairs (half,half)", list(itertools.combinations(half,2))),
                   ("pairs from cols+halfs", list(itertools.combinations(cols+half,2)))]:
    info=[]
    gf = family_gens(ring, A, fams, m, cache=cache, info=info)
    Jl = ideal_dims_modp(ring, merge(gw, gf), m)
    print(f"  J_w + {name} [{len(fams)} fams]: {Jl} cum_eq={cum_eq(Jl,JA)}", flush=True)
    print("     types:", sorted(set((sz, rho, tuple(sorted(gi.items()))) for F, sz, rho, gi in info))[:5])
