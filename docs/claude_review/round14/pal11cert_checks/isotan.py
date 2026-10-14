"""Z_5-isotypic computation, mod a prime p = 1 mod 5, of
  (i)  rank_p K_{<=8} (K = evaluation matrix of all reduced monomials of degree <= 8 on A), and per-degree ranks,
  (ii) the tangent space T_8 = E_8(ker E_{<8}) in isotypic coordinates (rref bases), saved for membership tests.
G = Z_5 = <g>, g: row 1+i -> row 1+(3i mod 11) (coordinate x=i -> 3i mod 11 of k').  For a character chi_j(g^t) = w^(j t)
(w a primitive 5th root of unity mod p) the block is B_j[a,k] = sum_t chi_j(g^t)^(-1) E(a, g^t k) on free orbit reps,
plus (j = 0 only) rows/columns of fixed monomials/points with entry E(a,k).  rank_p E = sum_j rank_p B_j.
usage: python3 isotan.py p out.pkl"""
import sys, time, pickle
import numpy as np
from flint import nmod_mat
from fplib import pal_patterns, box, reduced_monomials
p = int(sys.argv[1]); out = sys.argv[2]
assert (p - 1) % 5 == 0
# primitive 5th root of unity
for h in range(2, 1000):
    w = pow(h, (p - 1) // 5, p)
    if w != 1: break
M, pats = pal_patterns(11); A = box(M, pats, 2); d = 11
perm = [(3 * i) % 11 for i in range(11)]          # coordinate i -> perm[i]
def act_pt(k):
    out_ = [0] * d
    for i in range(d): out_[perm[i]] = k[i]
    return tuple(out_)
def orbits(items):
    S = set(items); seen = set(); reps = []; fixed = []
    for it in items:
        if it in seen: continue
        orb = [it]
        for t in range(4): orb.append(act_pt(orb[-1]))
        assert act_pt(orb[-1]) == it
        for o in orb: assert o in S
        seen.update(orb)
        if len(set(orb)) == 1: fixed.append(it)
        else:
            assert len(set(orb)) == 5
            reps.append(orb)          # orb[t] = g^t . rep
    return reps, fixed
t0 = time.time()
Pre, Pfix = orbits(A)
print('points: free orbits', len(Pre), 'fixed', len(Pfix), flush=True)
Aidx = {k: i for i, k in enumerate(A)}
Aarr = np.array(A, dtype=np.int64)
def monoeval(a):  # values of monomial a on all of A (vector)
    col = np.ones(len(A), dtype=np.int64)
    for x, ax in enumerate(a):
        if ax == 1: col = col * Aarr[:, x]
        elif ax == 2: col = col * (Aarr[:, x] ** 2)
    return col
chi = [[pow(w, (j * t) % 5, p) for t in range(5)] for j in range(5)]
chiinv = [[pow(w, (-j * t) % 5, p) for t in range(5)] for j in range(5)]
Pidx = np.array([[Aidx[o] for o in orb] for orb in Pre], dtype=np.int64)   # nfree x 5
Fidx = np.array([Aidx[k] for k in Pfix], dtype=np.int64)
def block_rows(monos):
    """returns per-character list of row arrays (python ints mod p) for the given set of monomials (one degree)."""
    Mre, Mfix = orbits(monos)
    rows = [[] for _ in range(5)]
    for orb in Mre:
        a = orb[0]
        v = monoeval(a)
        V = v[Pidx]  # nfree x 5 : E(a, g^t k)
        for j in range(5):
            coef = np.array(chiinv[j], dtype=object)
            r = (V.astype(object) @ coef) % p
            if j == 0:
                r = list(r) + [int(x) % p for x in v[Fidx]]
            rows[j].append([int(x) for x in r])
    for a in Mfix:
        v = monoeval(a)
        V = v[Pidx]
        r = [(5 * int(x)) % p for x in V[:, 0]] + [int(x) % p for x in v[Fidx]]   # E(a, sum_t g^t k) = 5 E(a,k)
        rows[0].append(r)
    return rows, Mre, Mfix
ncol = [len(Pre) + len(Pfix)] + [len(Pre)] * 4
res = {'p': p, 'w': w, 'perm': perm, 'ncol': ncol}
low = [[] for _ in range(5)]
ranks_by_deg = []
for e in range(9):
    monos = reduced_monomials(d, e)
    rows, Mre, Mfix = block_rows(monos)
    if e < 8:
        for j in range(5): low[j] += rows[j]
        rk = []
        for j in range(5):
            if low[j]:
                rk.append(nmod_mat(len(low[j]), ncol[j], [x for r in low[j] for x in r], p).rank())
            else: rk.append(0)
        ranks_by_deg.append(rk)
        print('deg <=', e, 'block ranks', rk, 'total', sum(rk), '(%.0fs)' % (time.time() - t0), flush=True)
    else:
        top = rows; res['Mre8'] = [orb for orb in Mre]; res['Mfix8'] = Mfix
res['ranks_low'] = ranks_by_deg
# full rank of K_{<=8} and tangent spaces
T = []
KB = []
fullrk = []
for j in range(5):
    Klow = nmod_mat(len(low[j]), ncol[j], [x for r in low[j] for x in r], p)
    Kall = nmod_mat(len(low[j]) + len(top[j]), ncol[j], [x for r in low[j] + top[j] for x in r], p)
    fullrk.append(Kall.rank())
    X, nul = Klow.nullspace()
    # columns 0..nul-1 of X span the kernel
    Kb = nmod_mat(ncol[j], nul, [int(X[i, c]) for i in range(ncol[j]) for c in range(nul)], p)
    E8 = nmod_mat(len(top[j]), ncol[j], [x for r in top[j] for x in r], p)
    KB.append([[int(Kb[i, c]) for i in range(ncol[j])] for c in range(nul)])
    Tj = (E8 * Kb).transpose()   # nul x (#deg-8 rows of block j)
    R, rk = Tj.rref()
    T.append([[int(R[i, c]) for c in range(R.ncols())] for i in range(rk)])
    print('block', j, 'cols', ncol[j], 'kernel dim (deg<8)', nul, 'dim T_8 block', rk, 'rank K_<=8 block', fullrk[-1],
          '(%.0fs)' % (time.time() - t0), flush=True)
res['T'] = T; res['KB'] = KB; res['Pre'] = Pre; res['Pfix'] = Pfix; res['fullrank'] = fullrk
print('rank_p K_{<=8} =', sum(fullrk), ' |A| =', len(A), ' dim T_8 =', sum(len(t) for t in T))
pickle.dump(res, open(out, 'wb'))
