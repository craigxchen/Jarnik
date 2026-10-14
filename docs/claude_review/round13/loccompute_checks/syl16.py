"""Syl(16, L) with the Singer cycle Z_15 (coordinates = alpha^l, l=0..14, F_16 = F_2[t]/(t^4+t+1))."""
import sys, time, itertools
import numpy as np
from loclib import PR, image, pullback
from bigring import BigRing
from cyclib import CycRing, ideal_dims_cyc, hilbert_cyc
from hilbert_lean import hilbert_lean
from jfast2 import jproj_gens_fast2
import models

def gf16_powers():
    x = 1; out = []
    for l in range(15):
        out.append(x)
        x <<= 1
        if x & 16: x ^= 0b10011
    return out

def syl16_singer(L):
    pts = gf16_powers()             # coordinate l <-> field element pts[l]
    cols = []; colA = []
    for a in range(1, 16):
        cols.append([1 if bin(a & x).count('1') % 2 == 1 else 0 for x in pts]); colA.append(a)
    return 15, cols, [L]*15, colA

def main(L, m, which):
    t0 = time.time()
    d, cols, Ls, colA = syl16_singer(L)
    A = models.box(d, cols, Ls)
    print(f"Syl(16,{L}) |A|={len(A)} [{time.time()-t0:.0f}s]", flush=True)
    ring = BigRing(d, m); cr = CycRing(ring, 15)
    # check column set invariant under shift
    cs = {tuple(c) for c in cols}
    assert all(tuple(np.roll(c, 1)) in cs for c in cols)
    H, pc = hilbert_lean(cr, A, m)
    JA = [ring.N[i]-H[i] for i in range(m+1)]
    print(f"  H_p={H}  dimJ(A) in R={JA} [{time.time()-t0:.0f}s]", flush=True)
    gens = {}
    def addpow(v, tag):
        s = len(set((A @ np.array(v)).tolist()))
        if s <= m:
            vec = ring.one()
            for t in range(s): vec = ring.mul_lin(vec, t, v)
            if vec.any(): gens.setdefault(s, []).append(vec)
        print(f"    width dir {tag}: |image|={s}", flush=True)
    addpow(cols[0], 'column')
    addpow([1]*15, 'all-ones (=-k_base)')
    lines = sorted({tuple(sorted((a, b, a ^ b))) for a in range(1, 16) for b in range(1, 16) if a != b})
    idx = {a: i for i, a in enumerate(colA)}
    # half-sums of lines
    for ln in lines:
        z = (np.array(cols[idx[ln[0]]]) + cols[idx[ln[1]]] + cols[idx[ln[2]]]) // 2
        addpow(list(z), f'half-sum line {ln}')
    for c in range(1, 16):
        V = [a for a in range(1, 16) if bin(a & c).count('1') % 2 == 0]
        z = sum(np.array(cols[idx[a]]) for a in V); assert (z % 4 == 0).all()
        addpow(list(z // 4), f'quarter-sum plane c={c}')
    for l in range(15):
        e = [0]*15; e[l] = 1; addpow(e, f'coordinate {l}')
    cum = lambda a: list(itertools.accumulate(a))
    def report(tag):
        G = {k: np.array(v) for k, v in gens.items()}
        Jl = ideal_dims_cyc(cr, G, m, verbose=True)
        print(f"  {tag}: {Jl}\n     cum_eq={[x==y for x,y in zip(cum(Jl),cum(JA))]} [{time.time()-t0:.0f}s]", flush=True)
    report("J_w(columns, all-ones, half-sums)")
    cache = {}
    if 'lines' in which:
        for ln in lines:      # all lines (orbit reduction not needed for correctness)
            F = [cols[idx[a]] for a in ln]
            P = image(A, F)
            jp, rho, fails = jproj_gens_fast2(P, 3, m)
            for i, forms in jp.items():
                for f in forms:
                    vec = pullback(ring, F, f, i, PR, cache)
                    if vec.any(): gens.setdefault(i, []).append(vec)
        print(f"    line image |P|={len(P)} rho={rho} newgens={[(i,len(f)) for i,f in jp.items() if f]}", flush=True)
        report("+ line triples")
    if 'planes' in which:
        # 3-dim subspaces V of F_2^4 (as a's): 15 of them; family = 7 columns a in V\0
        planes = []
        for c in range(1, 16):   # V = {a : a.c = 0}
            planes.append([a for a in range(1, 16) if bin(a & c).count('1') % 2 == 0])
        for V in planes[:15]:
            F = [cols[idx[a]] for a in V]
            P = image(A, F)
            jp, rho, fails = jproj_gens_fast2(P, 7, m)
            print(f"    plane {V}: |P|={len(P)} rho={rho} newgens={[(i,len(f)) for i,f in jp.items() if f]} fails={fails} [{time.time()-t0:.0f}s]", flush=True)
            for i, forms in jp.items():
                for f in forms:
                    vec = pullback(ring, F, f, i, PR, cache)
                    if vec.any(): gens.setdefault(i, []).append(vec)
        report("+ plane 7-col families")

if __name__ == '__main__':
    main(int(sys.argv[1]), int(sys.argv[2]), sys.argv[3])
