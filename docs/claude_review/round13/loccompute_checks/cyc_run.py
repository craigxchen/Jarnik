"""Paley-12 (Z_11) / Syl(16) Singer (Z_15) locality runs."""
import sys, time, itertools, math
import numpy as np
from loclib import Ring, PR, image, jproj_gens, pullback
from cyclib import CycRing, ideal_dims_cyc, hilbert_cyc
import models

def dir_orbit_reps(d, R=1):
    seen = set(); reps = []
    for v in itertools.product(range(-R, R+1), repeat=d):
        nz = [t for t in v if t]
        if not nz or math.gcd(*[abs(t) for t in nz]) != 1: continue
        if v in seen: continue
        orb = set()
        for sgn in (1, -1):
            for t in range(d):
                orb.add(tuple(sgn*z for z in np.roll(v, t)))
        seen |= orb; reps.append(v)
    return reps

def image_sizes(A, V, lo, hi):
    A = np.asarray(A, dtype=np.int32)
    out = np.zeros(len(V), dtype=np.int64)
    for s in range(0, len(V), 2000):
        P = A @ np.asarray(V[s:s+2000], dtype=np.int32).T
        cnt = np.zeros(P.shape[1], dtype=np.int64)
        for t in range(lo, hi+1):
            cnt += (P == t).any(axis=0)
        out[s:s+2000] = cnt
    return out

def powgen(ring, v, e):
    vec = ring.one()
    for t in range(e):
        vec = ring.mul_lin(vec, t, v)
    return vec

def fam_orbit_reps(cols, b, q):
    seen = set(); reps = []
    for S in itertools.combinations(range(len(cols)), b):
        if S in seen: continue
        orb = {tuple(sorted((x + t) % q for x in S)) for t in range(q)}
        seen |= orb; reps.append(S)
    return reps

def main(model, L, m, bmax, R=1):
    t0 = time.time()
    if model == 'paley':
        d, cols, Ls = models.paley12(L); q = 11
    A = models.box(d, cols, Ls)
    ring = Ring(d, m); cr = CycRing(ring, q)
    print(f"{model} L={L} m={m} d={d} |A|={len(A)}  [{time.time()-t0:.0f}s]", flush=True)
    H = hilbert_cyc(cr, A, m, None)
    JA = [ring.N[i] - H[i] for i in range(m+1)]
    print(f"  H_p={H} sum={sum(H)}  dimJ(A) in R={JA}  [{time.time()-t0:.0f}s]", flush=True)
    V = dir_orbit_reps(d, R)
    sz = image_sizes(A, V, -2*d, 2*d)
    gens = {}
    nw = 0
    for v, s in zip(V, sz):
        if s <= m:
            g = powgen(ring, v, int(s))
            if g.any():
                gens.setdefault(int(s), []).append(g); nw += 1
    print(f"  width generators (orbit reps, R={R}): {nw} by degree {sorted((k, len(x)) for k, x in gens.items())} [{time.time()-t0:.0f}s]", flush=True)
    cum = lambda a: list(itertools.accumulate(a))
    G = {k: np.array(v) for k, v in gens.items()}
    Jw = ideal_dims_cyc(cr, G, m, verbose=True)
    print(f"  b=1 J_w dims={Jw}\n        cum_eq={[x==y for x,y in zip(cum(Jw), cum(JA))]} [{time.time()-t0:.0f}s]", flush=True)
    cache = {}
    for b in range(2, bmax+1):
        types = {}
        for S in fam_orbit_reps(cols, b, q):
            F = [cols[x] for x in S]
            P = image(A, F)
            jp, rho = jproj_gens(P, b, m)
            key = (len(P), rho, tuple((i, len(f)) for i, f in sorted(jp.items()) if f))
            types[key] = types.get(key, 0) + 1
            for i, forms in jp.items():
                for f in forms:
                    vec = pullback(ring, F, f, i, PR, cache)
                    if vec.any(): gens.setdefault(i, []).append(vec)
        print(f"  {b}-col family orbit types (|P|,rho,newgens)->#orbits: {types} [{time.time()-t0:.0f}s]", flush=True)
        G = {k: np.array(v) for k, v in gens.items()}
        Jl = ideal_dims_cyc(cr, G, m, verbose=True)
        print(f"  b={b} dims={Jl}\n        cum_eq={[x==y for x,y in zip(cum(Jl), cum(JA))]} [{time.time()-t0:.0f}s]", flush=True)
        if cum(Jl)[-1] == cum(JA)[-1] and all(x==y for x,y in zip(cum(Jl), cum(JA))):
            break

if __name__ == '__main__':
    main(sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4]))
