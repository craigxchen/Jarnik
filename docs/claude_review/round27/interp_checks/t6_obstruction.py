"""Compute the obstruction space O_e = (sum_l M_e(A_l))^perp / M_e(A)^perp for a model box, and test candidate forms.
Pairing: form F = sum_alpha F_alpha k^alpha (deg e) with moment vector m_alpha(c) = sum_k c_k k^alpha."""
import sys, itertools, time
sys.path.insert(0, '/tmp/claude-0/-home-user-Jarnik/ba4929a4-8665-5790-985a-3f09130a14b9/scratchpad/round27/interp_checks')
from shiftlib import *
from models import *
from flint import nmod_mat
name = sys.argv[1]; Lh = int(sys.argv[2]); e = int(sys.argv[3])
if name == 'sdc7': M, pats = sdc_patterns(7)
elif name == 'syl8': M, pats = syl_patterns(3)
A, inds = box(M, pats, Lh)
d = M - 1
Aset = set(A)
L = []
for x in range(d):
    for s in (1, -1):
        l = [0]*d; l[x] = s; L.append(tuple(l))
for x in range(d):
    for y in range(x+1, d):
        for s in (1, -1):
            for t in (1, -1):
                l = [0]*d; l[x] = s; l[y] = t; L.append(tuple(l))
Ne = monos_deg(d, e)
MA = M_space(A, e, d)
acc = None
for l in L:
    Al = [k for k in A if tuple(ki + li for ki, li in zip(k, l)) in Aset]
    Ml = M_space(Al, e, d)
    if Ml is None: continue
    acc = stack([acc, Ml], len(Ne)) if acc is not None else Ml
    R = acc.rref()[0]; rr = acc.rank()
    acc = nmod_mat(rr, len(Ne), [int(R[i, j]) for i in range(rr) for j in range(len(Ne))], P)
print('dim M_e(A)', MA.rank(), 'dim sum', acc.rank())
# kernel of acc (forms killing sum) and kernel of MA (= J(A)_e)
Ka, na = acc.nullspace()
KA, nA = MA.nullspace()
print('dim (sum)^perp', na, 'dim J(A)_e', nA)
# candidate forms: products of column forms l_p(k) = sum_x ind_p[x] k_x  (k(P)), and other symmetric candidates
import numpy as np
def poly_to_vec(coeffs):  # dict alpha->coef
    idx = {a: i for i, a in enumerate(Ne)}
    v = [0]*len(Ne)
    for a, c in coeffs.items():
        v[idx[a]] = c % P
    return v
def mul(p, q):
    r = {}
    for a, x in p.items():
        for b, y in q.items():
            c = tuple(i + j for i, j in zip(a, b)); r[c] = r.get(c, 0) + x*y
    return r
def lin(vec):
    return {tuple(1 if i == x else 0 for i in range(d)): vec[x] for x in range(d) if vec[x]}
def in_space(v, K, nk):
    # is v in span of columns of K (first nk columns)?
    rows = [[int(K[i, j]) for j in range(nk)] + [v[i]] for i in range(K.nrows())]
    m1 = nmod_mat(K.nrows(), nk, [x for r in rows for x in r[:-1]], P)
    m2 = nmod_mat(K.nrows(), nk + 1, [x for r in rows for x in r], P)
    return m1.rank() == m2.rank()
cols = [lin(ind) for ind in inds]
results = []
for S in itertools.combinations(range(len(cols)), e):
    p = {tuple([0]*d): 1}
    for i in S: p = mul(p, cols[i])
    v = poly_to_vec(p)
    a = in_space(v, Ka, na); b = in_space(v, KA, nA)
    results.append((S, a, b))
good = [r for r in results if r[1] and not r[2]]
print('products of', e, 'distinct column forms: total', len(results), '; in cap J(A_l) but not J(A):', len(good))
print(good[:10])
