"""Validation of the isotypic engine and of memlib (mod p2):
 1. map isotypic kernel vectors back to functions c on A; check E_{<8} c = 0 directly (all 17469 reduced monomials
    of degree < 8, no symmetry used);
 2. Phi_c coefficient vector mu_8(c) computed directly must be a member (memlib), for several c in each block;
 3. random vectors and perturbed mu_8(c) must be non-members;
 4. the 436 vectors c (one per kernel basis element) are linearly independent mod p (dim V(A,8) >= 436 mod p)."""
import pickle, random
import numpy as np
from flint import nmod_mat
from fplib import pal_patterns, box, reduced_monomials
from memlib import T8
fn = 'T8_p2.pkl'
r = pickle.load(open(fn, 'rb')); p = r['p']; w = r['w']
M_, pats = pal_patterns(11); A = box(M_, pats, 2); d = 11
Aidx = {k: i for i, k in enumerate(A)}
Pre, Pfix = r['Pre'], r['Pfix']
chiinv = [[pow(w, (-j * t) % 5, p) for t in range(5)] for j in range(5)]
def lift(j, kv):  # kernel vector in block j -> function on A (python ints mod p)
    c = [0] * len(A)
    for o, orb in enumerate(Pre):
        for t in range(5):
            c[Aidx[orb[t]]] = (c[Aidx[orb[t]]] + kv[o] * chiinv[j][t]) % p
    if j == 0:
        for f, k in enumerate(Pfix):
            c[Aidx[k]] = (c[Aidx[k]] + kv[len(Pre) + f]) % p
    return c
Aarr = np.array(A, dtype=np.int64)
def monovals(a):
    col = np.ones(len(A), dtype=np.int64)
    for x, ax in enumerate(a):
        if ax == 1: col = col * Aarr[:, x]
        elif ax == 2: col = col * Aarr[:, x] ** 2
    return col
low = [a for e in range(8) for a in reduced_monomials(d, e)]
top = reduced_monomials(d, 8)
Elow = np.array([monovals(a) for a in low], dtype=np.int64)    # 17469 x 12093, entries -1,0,1
Etop = np.array([monovals(a) for a in top], dtype=np.int64)
def matvec(E, c):  # exact mod p using split into 30-bit limbs
    c = np.array(c, dtype=object)
    lo = np.array([x & ((1 << 30) - 1) for x in c], dtype=np.int64)
    hi = np.array([x >> 30 for x in c], dtype=np.int64)
    a = E @ lo; b = E @ hi   # |entries| <= 12093 * 2^30 < 2^44: fine
    return [(int(x) + (int(y) << 30)) % p for x, y in zip(a, b)]
T = T8(fn)
random.seed(1)
allc = []
for j in range(5):
    KB = r['KB'][j]
    for idx in range(len(KB)):
        c = lift(j, KB[idx]); allc.append(c)
    for idx in random.sample(range(len(KB)), 3):
        c = lift(j, KB[idx])
        z = matvec(Elow, c)
        assert not any(z), 'kernel check failed'
        mu = matvec(Etop, c)
        vec = {a: v for a, v in zip(top, mu) if v}
        ok = T.member(vec)
        vec2 = dict(vec); a0 = top[random.randrange(len(top))]; vec2[a0] = (vec2.get(a0, 0) + 1) % p
        print('block', j, 'kernel elt', idx, ': E_<8 c = 0 OK; mu_8(c) member:', ok, '; perturbed member:', T.member(vec2), flush=True)
rv = {a: random.randrange(p) for a in top}
print('random vector member:', T.member(rv))
# independence of the lifted c's
Mc = nmod_mat(len(allc), len(A), [x for c in allc for x in c], p)
print('rank of the', len(allc), 'lifted kernel vectors:', Mc.rank())
