"""Standard-monomial version for subsets of the ternary cube {-1,0,1}^d (k_x^3 = k_x there, so Poly_{<e}|_S is
spanned by monomials with exponents <= 2 of degree < e, and on V(S,e) the degree-e moments of non-standard
monomials vanish).  Own code."""
import itertools
import numpy as np, flint
from klib import shiftset
from klib2 import evalmat, nmod
_cache = {}
def smonos(d, e):
    key = (d, e)
    if key not in _cache:
        out = []
        for a in itertools.product(range(3), repeat=d):
            if sum(a) == e: out.append(a)
        _cache[key] = out
    return _cache[key]
def smonos_lt(d, e):
    out = []
    for t in range(e): out += smonos(d, t)
    return out
def Me_modp(S, e, p, d):
    if len(S) == 0: return None
    lt = smonos_lt(d, e)
    if not lt:
        nul = len(S)
        K = flint.nmod_mat(len(S), nul, [1 if i == j else 0 for i in range(len(S)) for j in range(nul)], p)
    else:
        M = nmod(evalmat(S, lt, p), p)
        K, nul = M.nullspace()
        if nul == 0: return None
        K = flint.nmod_mat(len(S), nul, [int(K[i, j]) for i in range(len(S)) for j in range(nul)], p)
    E = nmod(evalmat(S, smonos(d, e), p), p)
    return E * K
def shift_dims(S, L, e, p, d, verbose=False):
    MS = Me_modp(S, e, p, d)
    hs = 0 if MS is None else MS.rank()
    nr = len(smonos(d, e))
    acc = None; cur = 0
    for i, l in enumerate(L):
        Sl = shiftset(S, l)
        if len(Sl) <= e: continue
        Ml = Me_modp(Sl, e, p, d)
        if Ml is None: continue
        if acc is None: acc = Ml
        else:
            ent = []
            for r in range(nr):
                ent += [int(acc[r, j]) for j in range(acc.ncols())] + [int(Ml[r, j]) for j in range(Ml.ncols())]
            acc = flint.nmod_mat(nr, acc.ncols() + Ml.ncols(), ent, p)
        T = acc.transpose().rref()[0]; r = acc.rank()
        acc = flint.nmod_mat(nr, r, [int(T[j, i]) for i in range(nr) for j in range(r)], p)
        cur = r
        if verbose and i % 200 == 0: print("   shift", i, "of", len(L), "dim sum", cur, "of", hs, flush=True)
        if cur == hs: break
    return hs, cur
