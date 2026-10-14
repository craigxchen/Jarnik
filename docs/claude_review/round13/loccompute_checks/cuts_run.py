"""Two- and three-cut blocks: J(A) vs J_w (b=1) and + pullbacks of column pairs (b=2), all columns (b=3).
Reports, per block, rho, first degree where J_w fails, minimal b (columns) for equality in degrees <= rho (+1)."""
import sys, time, itertools, json
import numpy as np
from locrun import *

def analyse(d, cols, Ls, R, extra_deg=1):
    A = models.box(d, cols, Ls)
    if len(A) == 0: return None
    ring0 = Ring(d, 2*d)
    H = hilbert_modp(A, ring0, 2*d)
    rho = max(i for i, h in enumerate(H) if h > 0)
    m = min(rho + extra_deg, 2*d)
    ring = Ring(d, m)
    JA = [ring.N[i]-H[i] for i in range(m+1)]
    gw = width_gens(ring, width_dirs(A, d, R, m))
    res = {'|A|': len(A), 'rho': rho, 'H': H[:rho+1]}
    gall = gw
    for b in range(1, len(cols)+1):
        if b > 1:
            gall = merge(gall, family_gens(ring, A, list(itertools.combinations(cols, b)), m))
        Jl = ideal_dims_modp(ring, gall, m)
        eq = cum_eq(Jl, JA)
        fail = [i for i in range(m+1) if not eq[i]]
        res[f'b{b}_firstfail'] = fail[0] if fail else None
        if not fail:
            res['b0'] = b; break
        if fail[0] <= rho: continue
        # fails only above rho
        res.setdefault('b0_upto_rho', b)
    res.setdefault('b0', None)
    if 'b0_upto_rho' not in res and res['b0'] is not None: res['b0_upto_rho'] = res['b0']
    return res

if __name__ == '__main__':
    which = sys.argv[1]; dmax = int(sys.argv[2])
    import os
    fn = f'cuts_{which}_d{dmax}_out.txt'
    done = set()
    if os.path.exists(fn):
        done = {l[:l.index('{')].strip() for l in open(fn) if '{' in l}
    out = open(fn, 'a')
    t0 = time.time()
    if which in ('two', 'two_a1'):
        for a, b, c in itertools.product(range(0, dmax+1), repeat=3):
            if a + b + c > dmax or b == 0 or a > c or a + b < 2 or b + c < 2: continue
            if which == 'two_a1' and (a == 0 or a + b + c < 6): continue
            d = a + b + c
            S1 = list(range(0, a+b)); S2 = list(range(a, a+b+c))
            for L1 in range(1, len(S1)):
                for L2 in range(1, len(S2)):
                    if a == c and L1 > L2: continue
                    key = f"two (a,b,c)=({a},{b},{c}) L=({L1},{L2})"
                    if key in done: continue
                    r = analyse(d, *models.cuts(d, [S1, S2], [L1, L2])[1:], R=2 if d <= 6 else 1)
                    if r is None: continue
                    line = key + " " + json.dumps(r)
                    print(line, flush=True); out.write(line + "\n"); out.flush()
    if which == 'three':
        regions = [(1,0,0),(0,1,0),(0,0,1),(1,1,0),(1,0,1),(0,1,1),(1,1,1)]
        seen = set()
        for sz in itertools.product(range(0, dmax+1), repeat=7):
            if sum(sz) > dmax: continue
            # canonical under S_3
            def perm(sz, p):
                d = {tuple(r[p[i]] for i in range(3)): z for r, z in zip(regions, sz)}
                return tuple(d[r] for r in regions)
            canon = min(perm(sz, p) for p in itertools.permutations(range(3)))
            if canon in seen: continue
            seen.add(canon)
            sz = canon
            sets = [[], [], []]; x = 0
            for r, z in zip(regions, sz):
                for _ in range(z):
                    for j in range(3):
                        if r[j]: sets[j].append(x)
                    x += 1
            d = x
            if any(len(S) < 2 for S in sets) or d < 3: continue
            if len({tuple(S) for S in sets}) < 3: continue
            for Lv in itertools.product((1, 2), repeat=3):
                if any(Lv[j] >= len(sets[j]) for j in range(3)): continue
                r = analyse(d, *models.cuts(d, sets, Lv)[1:], R=2 if d <= 6 else 1)
                if r is None: continue
                line = f"three regions={sz} L={Lv} " + json.dumps(r)
                out.write(line + "\n"); out.flush()
                if r.get('b0') != 1: print(line, flush=True)
    print(f"done [{time.time()-t0:.0f}s]")
