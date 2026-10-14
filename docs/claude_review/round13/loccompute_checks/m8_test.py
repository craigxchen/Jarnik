import sys, time, itertools
from locrun import *
k, L, m = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
d, cols, Ls = models.sylvester(k, L)
A = models.box(d, cols, Ls)
ring = Ring(d, m)
t0=time.time()
H = hilbert_modp(A, ring, m)
JA = [ring.N[i]-H[i] for i in range(m+1)]
print(f"Syl({2**k},{L}) |A|={len(A)} m={m} H_p={H} dimJ(A)_p={JA}", flush=True)
for R in (1,2):
    wd = width_dirs(A, d, R, m)
    gw = width_gens(ring, wd)
    Jw = ideal_dims_modp(ring, gw, m)
    print(f"  J_w(R={R}) #dirs={len(wd)} dims={Jw} cum_eq={cum_eq(Jw,JA)}  [{time.time()-t0:.0f}s]", flush=True)
cache={}
for b in (1,2,3):
    info=[]
    gf = family_gens(ring, A, list(itertools.combinations(cols, b)), m, cache=cache, info=info)
    Jl = ideal_dims_modp(ring, merge(gw, gf), m)
    print(f"  J_w + {b}-col pullbacks: dims={Jl} cum_eq={cum_eq(Jl,JA)}  [{time.time()-t0:.0f}s]", flush=True)
    rhos = sorted(set((len(F), sz, rho, tuple(sorted(gi.items()))) for F, sz, rho, gi in info))
    for r in rhos[:6]: print('     fam: b,|P|,rho,newgens', r)
