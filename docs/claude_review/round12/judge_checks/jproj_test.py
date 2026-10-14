"""J6b: add to J_w the pullbacks of the leading-form ideals J(pi_F(A)) of 2-dim projections (pairs of directions F)
and test J(A) = J_w + J_proj2 in degrees <= m (same mod-p sandwich; 2D nullspaces computed exactly over Q)."""
import itertools, math, sys, time
import numpy as np
from flint import nmod_mat, fmpz_mat, fmpz
sys.path.insert(0, '.')
from jw_test import PR, eval_rank_profile, widths, mono_tuple, sylvester8

def lin(v, d):
    return {tuple(1 if y == x else 0 for y in range(d)): v[x] % PR for x in range(d) if v[x]}

def mul(p, q):
    r = {}
    for a, ca in p.items():
        for b, cb in q.items():
            k = tuple(s+t for s, t in zip(a, b)); r[k] = (r.get(k, 0) + ca*cb) % PR
    return {k: c for k, c in r.items() if c}

def powl(v, e, d, cache):
    key = (v, e)
    if key not in cache:
        cache[key] = {tuple([0]*d): 1} if e == 0 else mul(powl(v, e-1, d, cache), lin(v, d))
    return cache[key]

def proj_leading_forms(A, v1, v2, D):
    """exact: top-degree parts (degree i) of I(Lambda)_{<=i}, Lambda = image of A under (v1,v2), i <= D.
    returns dict i -> list of dicts {(a,b): fmpq} (binary forms of degree i) spanning J(Lambda)_i"""
    Av = np.array(A, dtype=np.int64)
    L = sorted(set(zip((Av @ np.array(v1)).tolist(), (Av @ np.array(v2)).tolist())))
    out = {}
    for i in range(1, D+1):
        mons = [(a, t-a) for t in range(i+1) for a in range(t+1)]
        M = fmpz_mat(len(L), len(mons), [p[0]**a * p[1]**b for p in L for (a, b) in mons])
        X, nul = M.nullspace()
        top = [(a, b) for (a, b) in mons if a + b == i]
        forms = []
        for c in range(nul):
            f = {}
            for r, (a, b) in enumerate(mons):
                if a + b == i and X[r, c] != 0: f[(a, b)] = X[r, c]
            if f: forms.append(f)
        out[i] = (forms, len(L))
    return out

def jdims(gens_by_deg, d, m):
    dims = []; prev = []
    for i in range(m+1):
        mi = [mono_tuple(c, d) for c in itertools.combinations_with_replacement(range(d), i)]
        idx = {e: t for t, e in enumerate(mi)}
        rows = []
        for r in prev:
            for x in range(d):
                rows.append({tuple(e[y] + (1 if y == x else 0) for y in range(d)): c for e, c in r.items()})
        rows += gens_by_deg.get(i, [])
        if not rows: dims.append(0); prev = []; continue
        M = nmod_mat(len(rows), len(mi), PR)
        for a, r in enumerate(rows):
            for e, c in r.items(): M[a, idx[e]] = c
        Rr, rk = M.rref()
        dims.append(rk)
        prev = [{mi[b]: int(Rr[a, b]) for b in range(len(mi)) if int(Rr[a, b])} for a in range(rk)]
    return dims

def modp(q):
    return int(q) % PR

def run(name, A, d, m, pairdirs, R=3):
    t0 = time.time()
    ranks = eval_rank_profile(A, d, m)
    H = [ranks[0]] + [ranks[i]-ranks[i-1] for i in range(1, m+1)]
    N = [math.comb(i+d-1, d-1) for i in range(m+1)]
    JA = [N[i]-H[i] for i in range(m+1)]
    wg = widths(A, d, R, m-1)
    cache = {}
    g = {}
    for v, h in wg:
        g.setdefault(h+1, []).append(powl(v, h+1, d, cache))
    Jw = jdims(g, d, m)
    # pairs
    g2 = {k: list(vs) for k, vs in g.items()}
    npairs = 0
    for v1, v2 in itertools.combinations(pairdirs, 2):
        npairs += 1
        lf = proj_leading_forms(A, v1, v2, m)
        for i, (forms, _) in lf.items():
            for f in forms:
                poly = {}
                for (a, b), q in f.items():
                    t = mul(powl(v1, a, d, cache), powl(v2, b, d, cache))
                    cq = modp(q)
                    for e, c in t.items(): poly[e] = (poly.get(e, 0) + c*cq) % PR
                poly = {e: c for e, c in poly.items() if c}
                if poly: g2.setdefault(i, []).append(poly)
    Jp = jdims(g2, d, m)
    cA = list(itertools.accumulate(JA))
    s = (f"{name}: |A|={len(A)} m={m} width-gens(R={R})={len(wg)} pairs={npairs}\n  dimJ(A)      ={JA}\n  dimJ_w       ={Jw}\n"
         f"  dimJ_w+proj2 ={Jp}\n  J(A)=J_w (cum)       : {[a==b for a,b in zip(cA, itertools.accumulate(Jw))]}\n"
         f"  J(A)=J_w+proj2 (cum) : {[a==b for a,b in zip(cA, itertools.accumulate(Jp))]}   [{time.time()-t0:.0f}s]")
    print(s, flush=True); return s

if __name__ == '__main__':
    which = sys.argv[1]
    out = []
    if which == 'M8':
        A = sylvester8(2)
        def dot(a, x): return bin(a & x).count('1') % 2
        others = list(range(7)); base = 7
        cols = [tuple((1 if dot(a, x) == 0 else 0) - (1 if dot(a, base) == 0 else 0) for x in others) for a in range(1, 8)]
        out.append(run('M8 pairs of columns', A, 7, 8, cols, R=3))
    open(f"jproj_test_{which}_out.txt", 'w').write("\n".join(out) + "\n")
