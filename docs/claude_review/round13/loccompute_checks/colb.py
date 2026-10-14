"""For a model: J(A) vs J_w(R) + pullbacks of all <=b column families, b = 1.. bmax (mod-p sandwich)."""
import sys, time, itertools, json
import numpy as np
from locrun import *

def run(name, d, cols, Ls, m, R=2, bmax=None, extra_dirs=(), out=None):
    t0 = time.time()
    A = models.box(d, cols, Ls)
    ring = Ring(d, m)
    H = hilbert_modp(A, ring, m); JA = [ring.N[i]-H[i] for i in range(m+1)]
    s = [f"{name}: d={d} |A|={len(A)} m={m} H_p={H}  dimJ(A) (in R)={JA}"]
    wd = width_dirs(A, d, R, m)
    gw = width_gens(ring, wd)
    Jw = ideal_dims_modp(ring, gw, m)
    s.append(f"  b=1 (J_w, R={R}, {len(wd)} dirs): {Jw} cum_eq={cum_eq(Jw,JA)}")
    print("\n".join(s), flush=True)
    res = {1: Jw}
    cache = {}
    gall = dict(gw)
    bmax = bmax or len(cols)
    for b in range(2, bmax+1):
        info = []
        gf = family_gens(ring, A, list(itertools.combinations(cols, b)), m, cache=cache, info=info)
        gall = merge(gall, gf)
        Jl = ideal_dims_modp(ring, gall, m)
        res[b] = Jl
        types = {}
        for F, sz, rho, gi in info:
            key = (sz, rho, tuple(sorted(gi.items())))
            types[key] = types.get(key, 0) + 1
        line = f"  b={b} (+{b}-col pullbacks): {Jl} cum_eq={cum_eq(Jl,JA)}  types(|P|,rho,newgens)->count: {types} [{time.time()-t0:.0f}s]"
        print(line, flush=True); s.append(line)
        if all(cum_eq(Jl, JA)):
            break
    if out:
        open(out, 'a').write("\n".join(s) + "\n")
    return JA, res

if __name__ == '__main__':
    which = sys.argv[1]
    if which == 'syl8':
        for L, m in ((1, 3), (1, 4), (2, 8), (3, 11)):
            run(f"Syl(8,{L})", *models.sylvester(3, L), m, out='colb_syl8_out.txt')
