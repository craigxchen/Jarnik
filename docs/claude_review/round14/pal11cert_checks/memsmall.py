"""Membership of count-feasible products in J(A)_m^perp for the small round-13 models (d = 7), no symmetry.
usage: python3 memsmall.py MODEL L m prefix [maxleaves]
* rank_Q K_{<=m} (exact, fmpz_mat) and rank_p K_{<=m} (p = 2^31 - 1 - ...); negatives certified only if equal.
* T_m mod p = E_m(ker_p E_{<m}); membership of (P_a / multinom(m;a)) by rank.
* members: exact rational check (fmpq rank of [K | b] == rank K) and zonotope test (translate inside A)."""
import sys, itertools
from math import factorial
import numpy as np
from flint import nmod_mat, fmpz_mat, fmpq_mat, fmpq
from fplib import MODELS, box, reduced_monomials
from memlib import product_coeffs
MODEL, L, m, pre = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), sys.argv[4]
maxleaves = int(sys.argv[5]) if len(sys.argv) > 5 else 10 ** 9
p = 2147483629
M, pats = MODELS[MODEL](); A = box(M, pats, L); d = M - 1
Aarr = np.array(A, dtype=np.int64)
def mv(a):
    col = np.ones(len(A), dtype=np.int64)
    for x, ax in enumerate(a):
        if ax == 1: col = col * Aarr[:, x]
        elif ax == 2: col = col * Aarr[:, x] ** 2
    return col
low = [a for e in range(m) for a in reduced_monomials(d, e)]
top = reduced_monomials(d, m)
El = np.array([mv(a) for a in low], dtype=np.int64)
Et = np.array([mv(a) for a in top], dtype=np.int64)
K = np.vstack([El, Et])
rQ = fmpz_mat(K.shape[0], K.shape[1], [int(x) for x in K.flatten()]).rank()
rp = nmod_mat(K.shape[0], K.shape[1], [int(x) % p for x in K.flatten()], p).rank()
print(MODEL, 'L', L, 'm', m, '|A|', len(A), ' rank_Q K_{<=m} =', rQ, ' rank_p =', rp, ' (negatives rigorous:', rQ == rp, ')', flush=True)
Nl = nmod_mat(El.shape[0], El.shape[1], [int(x) % p for x in El.flatten()], p)
X, nul = Nl.nullspace()
Kb = nmod_mat(len(A), nul, [int(X[i, j]) for i in range(len(A)) for j in range(nul)], p)
T = (nmod_mat(Et.shape[0], Et.shape[1], [int(x) % p for x in Et.flatten()], p) * Kb).transpose()
Tr, rk = T.rref()
Trows = [[int(Tr[i, j]) for j in range(Tr.ncols())] for i in range(rk)]
print('dim T_m mod p =', rk, flush=True)
idx = {a: i for i, a in enumerate(top)}
def multinom(a):
    v = factorial(m)
    for ax in a: v //= factorial(ax)
    return v
pool = []
with open(pre + '_pool.txt') as f:
    f.readline(); f.readline()
    for line in f: pool.append([int(x) for x in line.split()[1:1 + M]])
Aset = set(A)
def zonotope(G):
    G = np.array(G)
    verts = np.array([G[list(I)].sum(axis=0) if I else np.zeros(d, dtype=np.int64)
                      for r in range(len(G) + 1) for I in itertools.combinations(range(len(G)), r)])
    lo = verts.min(axis=0); hi = verts.max(axis=0)
    for k0 in itertools.product(*[range(-1 - int(l), 2 - int(h)) for l, h in zip(lo, hi)]):
        if all(tuple(v + np.array(k0)) in Aset for v in verts): return True
    return False
nmem = nnon = nzon = nexact = 0
members = []
for li, line in enumerate(open(pre + '_leaves.txt')):
    if li >= maxleaves: break
    Z = [pool[int(i)][1:] for i in line.split()]
    poly = product_coeffs(Z)
    if any(max(e) > 2 for e in poly):
        nnon += 1; continue
    vec = [0] * len(top)
    for e, c in poly.items(): vec[idx[e]] = (c * pow(multinom(e), -1, p)) % p
    Mm = nmod_mat(rk + 1, len(top), [x for r in Trows for x in r] + vec, p)
    if Mm.rank() == rk:
        nmem += 1; members.append(Z)
    else: nnon += 1
print('leaves tested', nmem + nnon, ': members mod p', nmem, ', non-members', nnon, flush=True)
for Z in members:
    zon = zonotope(Z)
    nzon += zon
    if nexact < 5:
        poly = product_coeffs(Z)
        rhs = [fmpq(0)] * len(low) + [fmpq(0)] * len(top)
        for e, c in poly.items(): rhs[len(low) + idx[e]] = fmpq(c, multinom(e))
        Aug = fmpq_mat(K.shape[0], K.shape[1] + 1, [x for r in range(K.shape[0]) for x in [int(y) for y in K[r]] + [rhs[r]]])
        r2 = Aug.rref()[1]
        print('   member', Z, ' zonotope:', zon, ' exact rank [K|b] =', r2, ' == rank K:', r2 == rQ, flush=True)
        nexact += 1
print('members that are zonotopes:', nzon, 'of', len(members))
