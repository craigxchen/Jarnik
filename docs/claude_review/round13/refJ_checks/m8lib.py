"""Independent library for the M8 locality checks (referee of J6(c)).

M8 = Syl(8), L = 2, m = 8.  Points x in F_2^3 = {0..7}, base point 0, non-base points 1..7 (coordinate i <-> x = i+1).
k' in {-1,0,1}^7, k_base = -sum k'.  Sign column of a != 0: <k, beta_a> = sum_x k_x (-1)^{a.x} = 2 k(S_a),
S_a = {x : a.x = 0} (uses sum k = 0).  Since 0 in S_a, k(S_a) = -k'(T_a), T_a = {x != 0 : a.x = 1} (4 points).
Box: |<k,beta_a>| <= 2L  <=>  |k'(T_a)| <= L for all 7 a.  Column direction (in k'-space): c_a = 1_{T_a}
(sign irrelevant for ideals).
"""
import itertools, math
import numpy as np
from flint import fmpz_mat, nmod_mat

PR = 2**31 - 1
d = 7

def dot(a, x):
    return bin(a & x).count('1') & 1

def m8_box(L=2):
    T = {a: [i for i in range(7) if dot(a, i + 1) == 1] for a in range(1, 8)}
    A = [k for k in itertools.product((-1, 0, 1), repeat=7) if all(abs(sum(k[i] for i in T[a])) <= L for a in T)]
    cols = {a: tuple(1 if i in T[a] else 0 for i in range(7)) for a in T}
    return A, cols

_mon_cache = {}
def monos(nv, i):
    key = (nv, i)
    if key not in _mon_cache:
        out = []
        for c in itertools.combinations_with_replacement(range(nv), i):
            e = [0] * nv
            for x in c: e[x] += 1
            out.append(tuple(e))
        _mon_cache[key] = out
    return _mon_cache[key]

def mono_index(nv, i):
    return {e: t for t, e in enumerate(monos(nv, i))}

def pmul(p, q):
    r = {}
    for a, ca in p.items():
        for b, cb in q.items():
            k = tuple(x + y for x, y in zip(a, b))
            r[k] = r.get(k, 0) + ca * cb
    return {k: c for k, c in r.items() if c}

def lin(v):
    out = {}
    for i, c in enumerate(v):
        if c:
            e = [0] * len(v); e[i] = 1; out[tuple(e)] = c
    return out

def ppow(p, e, nv):
    r = {tuple([0] * nv): 1}
    for _ in range(e): r = pmul(r, p)
    return r

def eval_ranks(A, m, exact=True):
    """cumulative rank of evaluation of monomials of degree <= i (exponents <= 2 suffice on {-1,0,1}^7), i = 0..m"""
    Aa = np.array(A, dtype=np.int64)
    cols = []; ranks = []
    for i in range(m + 1):
        for e in monos(d, i):
            if max(e) <= 2:
                col = np.ones(len(A), dtype=np.int64)
                for x, ex in enumerate(e):
                    if ex: col = col * Aa[:, x] ** ex
                cols.append(col)
        M = np.array(cols, dtype=np.int64)  # rows = monomials
        if exact:
            ranks.append(fmpz_mat(M.shape[0], M.shape[1], [int(t) for t in M.flatten()]).rank())
        else:
            ranks.append(nmod_mat(M.shape[0], M.shape[1], [int(t) % PR for t in M.flatten()], PR).rank())
    return ranks

def hilbert_JA(A, m, exact=True):
    r = eval_ranks(A, m, exact)
    H = [r[0]] + [r[i] - r[i - 1] for i in range(1, m + 1)]
    return H, [math.comb(i + d - 1, d - 1) - H[i] for i in range(m + 1)]

def width_set(A, R=3, maxdeg=8):
    """primitive v in [-R,R]^7 (first nonzero > 0) with tau+1 <= maxdeg; returns list (v, tau, h)."""
    Aa = np.array(A, dtype=np.int64)
    vs = [v for v in itertools.product(range(-R, R + 1), repeat=7)
          if any(v) and [t for t in v if t][0] > 0 and math.gcd(*[abs(t) for t in v if t]) == 1]
    V = np.array(vs, dtype=np.int64)
    out = []
    for s in range(0, len(V), 50000):
        blk = V[s:s + 50000]
        P = Aa @ blk.T
        h = P.max(axis=0) - P.min(axis=0)
        for t in np.nonzero(h <= 3 * maxdeg)[0]:
            vals = np.unique(P[:, t])
            tau = len(vals) - 1
            if tau + 1 <= maxdeg:
                out.append((tuple(int(z) for z in blk[t]), tau, int(h[t])))
    return out

def proj_ideal_gens(A, F, maxdeg):
    """For the image Y = pi_F(A) in Z^q, return {e: list of integer forms (dict) in q vars} = a basis of
    complements of (J(Y)_{e-1} * y) in J(Y)_e, e = 1..maxdeg (exact over Q). Also rho(Y) and |Y|."""
    Aa = np.array(A, dtype=np.int64)
    Y = sorted(set(map(tuple, (Aa @ np.array(F, dtype=np.int64).T).tolist())))
    q = len(F)
    gens = {}; prev = []   # prev: basis of J(Y)_{e-1} as dicts
    rho = None
    allm = []
    for e in range(0, maxdeg + 1):
        allm = allm + monos(q, e)
        M = fmpz_mat(len(Y), len(allm), [int(np.prod([p[z] ** ex[z] for z in range(q)])) for p in Y for ex in allm])
        X, nul = M.nullspace()
        topm = monos(q, e); tix = mono_index(q, e)
        forms = []
        for c in range(nul):
            f = {ex: int(X[r, c]) for r, ex in enumerate(allm) if sum(ex) == e and X[r, c] != 0}
            if f: forms.append(f)
        # rank of J(Y)_e
        def rows_of(fs):
            return [[f.get(ex, 0) for ex in topm] for f in fs]
        JYe = forms
        rk = fmpz_mat(len(JYe), len(topm), [x for r in rows_of(JYe) for x in r]).rank() if JYe else 0
        if rk == len(topm) and rho is None and e > 0:
            rho = e - 1
        # products prev * y_j
        prods = []
        for f in prev:
            for j in range(q):
                ej = [0] * q; ej[j] = 1
                prods.append(pmul(f, {tuple(ej): 1}))
        base_rows = rows_of(prods)
        cur = base_rows[:]
        r0 = fmpz_mat(len(cur), len(topm), [x for r in cur for x in r]).rank() if cur else 0
        new = []
        for f in JYe:
            trial = cur + rows_of([f])
            r1 = fmpz_mat(len(trial), len(topm), [x for r in trial for x in r]).rank()
            if r1 > r0:
                cur = trial; r0 = r1; new.append(f)
        if e > 0 and new: gens[e] = new
        # basis of J(Y)_e for next step
        prev = prods_basis = None
        # choose a basis of J(Y)_e among JYe
        bas = []; cur2 = []; r2 = 0
        for f in JYe:
            trial = cur2 + rows_of([f])
            r3 = fmpz_mat(len(trial), len(topm), [x for r in trial for x in r]).rank()
            if r3 > r2: cur2 = trial; r2 = r3; bas.append(f)
        prev = bas
    return gens, rho, len(Y)

def pullback(f, F, cache):
    """substitute y_z = <x, F[z]> in the q-variate form f"""
    out = {}
    for ex, c in f.items():
        t = {tuple([0] * d): c}
        for z, ez in enumerate(ex):
            if ez:
                key = (F[z], ez)
                if key not in cache: cache[key] = ppow(lin(F[z]), ez, d)
                t = pmul(t, cache[key])
        for k, cc in t.items(): out[k] = out.get(k, 0) + cc
    return {k: c for k, c in out.items() if c}

def ideal_degree_rows(gens_by_deg, i):
    """rows (dicts) spanning (ideal)_i: g * monomial for each generator g of degree e <= i"""
    rows = []
    for e, gl in gens_by_deg.items():
        if e > i: continue
        ms = monos(d, i - e)
        for g in gl:
            for mm in ms:
                rows.append({tuple(a + b for a, b in zip(k, mm)): c for k, c in g.items()})
    return rows

def dim_exact(rows, i):
    ix = mono_index(d, i); N = len(ix)
    if not rows: return 0
    data = [0] * (len(rows) * N)
    for r, f in enumerate(rows):
        for k, c in f.items(): data[r * N + ix[k]] = c
    return fmpz_mat(len(rows), N, data).rank()

def dim_modp(rows, i):
    ix = mono_index(d, i); N = len(ix)
    if not rows: return 0
    data = [0] * (len(rows) * N)
    for r, f in enumerate(rows):
        for k, c in f.items(): data[r * N + ix[k]] = c % PR
    return nmod_mat(len(rows), N, data, PR).rank()
