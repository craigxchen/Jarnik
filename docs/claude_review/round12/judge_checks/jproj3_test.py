"""J6c: J(M8) versus J_w + pullbacks of J(pi_F(A)) for triples F of columns (and of other small directions)."""
import itertools, math, sys, time
import numpy as np
from flint import nmod_mat, fmpz_mat
sys.path.insert(0, '.')
from jw_test import PR, eval_rank_profile, widths, mono_tuple, sylvester8
from jproj_test import mul, powl, jdims

def proj_leading_forms_k(A, F, D):
    Av = np.array(A, dtype=np.int64)
    img = sorted(set(map(tuple, (Av @ np.array(F).T).tolist())))
    q = len(F); out = {}
    for i in range(1, D+1):
        mons = [e for t in range(i+1) for e in (mono_tuple(c, q) for c in itertools.combinations_with_replacement(range(q), t))]
        M = fmpz_mat(len(img), len(mons), [int(np.prod([p[z]**e[z] for z in range(q)])) for p in img for e in mons])
        X, nul = M.nullspace()
        forms = []
        for c in range(nul):
            f = {e: int(X[r, c]) for r, e in enumerate(mons) if sum(e) == i and X[r, c] != 0}
            if f: forms.append(f)
        out[i] = forms
    return out, len(img)

def pullback(f, F, d, cache):
    poly = {}
    for e, cq in f.items():
        t = {tuple([0]*d): 1}
        for z, ez in enumerate(e):
            if ez: t = mul(t, powl(F[z], ez, d, cache))
        for k, c in t.items(): poly[k] = (poly.get(k, 0) + c*(cq % PR)) % PR
    return {k: c for k, c in poly.items() if c}

if __name__ == '__main__':
    A = sylvester8(2); d = 7; m = 8
    def dot(a, x): return bin(a & x).count('1') % 2
    others = list(range(7)); base = 7
    cols = [tuple((1 if dot(a, x) == 0 else 0) - (1 if dot(a, base) == 0 else 0) for x in others) for a in range(1, 8)]
    ranks = eval_rank_profile(A, d, m)
    H = [ranks[0]] + [ranks[i]-ranks[i-1] for i in range(1, m+1)]
    JA = [math.comb(i+d-1, d-1) - H[i] for i in range(m+1)]
    cache = {}; g = {}
    for v, h in widths(A, d, 2, m-1):
        g.setdefault(h+1, []).append(powl(v, h+1, d, cache))
    out = []
    for q in (3, 4):
        g2 = {k: list(v) for k, v in g.items()}
        t0 = time.time(); nf = 0
        for F in itertools.combinations(cols, q):
            lf, sz = proj_leading_forms_k(A, F, m)
            for i, forms in lf.items():
                for f in forms:
                    p = pullback(f, F, d, cache)
                    if p: g2.setdefault(i, []).append(p); nf += 1
        Jp = jdims(g2, d, m)
        s = (f"M8: J_w + {q}-column projections ({math.comb(7,q)} families, {nf} forms): dims={Jp}\n   J(A) dims={JA}\n"
             f"   cumulative equal: {[a==b for a,b in zip(itertools.accumulate(JA), itertools.accumulate(Jp))]} [{time.time()-t0:.0f}s]")
        print(s, flush=True); out.append(s)
    open('jproj3_test_M8_out.txt', 'w').write("\n".join(out) + "\n")
