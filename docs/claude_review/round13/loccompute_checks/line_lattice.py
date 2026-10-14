"""Abstract line image P_L = {y in [-L,L]^3 : y1+y2+y3 even} (the Sylvester line family when every coset is large).
Compare J(P_L) with J_w(P_L) (width powers of all primitive directions of the dual lattice, coefficients <= R),
in degrees <= min(8L, rho+1).  Full polynomial ring in 3 variables (no x^3 quotient), mod-p sandwich."""
import itertools, math, sys, time
import numpy as np
from flint import nmod_mat
PR = 2147475331
def monos(b, i): return [tuple(sum(1 for z in c if z == q) for q in range(b)) for c in itertools.combinations_with_replacement(range(b), i)]
def hilb(P, D):
    P = np.array(P, dtype=object); ranks = []; rows = []
    for i in range(D + 1):
        for e in monos(3, i): rows.append([int(p[0]**e[0] * p[1]**e[1] * p[2]**e[2]) % PR for p in P])
        ranks.append(nmod_mat(len(rows), len(P), [x for r in rows for x in r], PR).rank())
    return [ranks[0]] + [ranks[i] - ranks[i-1] for i in range(1, D + 1)]
def powform(v, e):
    # (v1 y1 + v2 y2 + v3 y3)^e as dict over monos
    out = {}
    for a in range(e + 1):
        for b in range(e - a + 1):
            c = e - a - b
            coef = math.comb(e, a) * math.comb(e - a, b) * v[0]**a * v[1]**b * v[2]**c
            if coef: out[(a, b, c)] = coef
    return out
def ideal_dims(gens, D):
    dims = []; basis = []
    for i in range(D + 1):
        mi = monos(3, i); idx = {e: t for t, e in enumerate(mi)}
        rows = []
        for r in basis:
            for q in range(3):
                w = [0] * len(mi)
                for e, c in r.items():
                    e2 = list(e); e2[q] += 1; w[idx[tuple(e2)]] = c
                rows.append(w)
        for g in gens.get(i, []):
            w = [0] * len(mi)
            for e, c in g.items(): w[idx[e]] = c % PR
            rows.append(w)
        if not rows: dims.append(0); basis = []; continue
        M = nmod_mat(len(rows), len(mi), [x % PR for r in rows for x in r], PR)
        Rr, rk = M.rref(); dims.append(rk)
        ent = [int(x) for x in Rr.entries()]
        basis = [{mi[c]: ent[a * len(mi) + c] for c in range(len(mi)) if ent[a * len(mi) + c]} for a in range(rk)]
    return dims
R = int(sys.argv[1]) if len(sys.argv) > 1 else 3
for L in range(1, 8):
    P = [y for y in itertools.product(range(-L, L + 1), repeat=3) if sum(y) % 2 == 0]
    D = min(8 * L, 3 * 2 * L + 2)
    H = hilb(P, D)
    rho = max(i for i, h in enumerate(H) if h > 0)
    D = min(D, rho + 1)
    JP = [len(monos(3, i)) - H[i] for i in range(D + 1)]
    # dual lattice directions: integer w with all entries of equal parity, form = (w.y)/2 if all odd else w.y (primitive)
    gens = {}; nd = 0
    Pa = np.array(P)
    for w in itertools.product(range(-2 * R, 2 * R + 1), repeat=3):
        if not any(w) or len({x % 2 for x in w}) > 1: continue
        if all(x % 2 == 0 for x in w):
            v = tuple(x // 2 for x in w)
        else:
            v = w   # form (w.y)/2 ; use w with scaling (same pullback ideal up to scalar)
        g = math.gcd(*[abs(x) for x in w])
        if g != 1 and not all(x % 2 == 0 for x in w): continue
        if all(x % 2 == 0 for x in w) and math.gcd(*[abs(x) for x in v]) != 1: continue
        vals = set((Pa @ np.array(v)).tolist())
        s = len(vals)
        if s <= D:
            gens.setdefault(s, []).append(powform(v, s)); nd += 1
    Jw = ideal_dims(gens, D)
    cum = lambda a: list(itertools.accumulate(a))
    eq = [x == y for x, y in zip(cum(Jw), cum(JP))]
    fail = [i for i, e in enumerate(eq) if not e]
    print(f"L={L} |P|={len(P)} rho={rho} D={D} #dirs={nd}  J(P)={JP}\n      J_w={Jw}  first fail: {fail[0] if fail else None}  (8L={8*L})", flush=True)
