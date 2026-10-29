"""Fast mod-p least-space routines (numpy evaluation + flint nmod_mat). Own code."""
import itertools
import numpy as np
import flint
from klib import monos, monos_lt, shiftset, kfull, kS, build_box, Lpair, Lcoord

def evalmat(S, mons, p):
    P = np.array(S, dtype=np.int64)            # |S| x d
    Ex = np.array(mons, dtype=np.int64)        # N x d
    if len(mons) == 0:
        return np.zeros((0, len(S)), dtype=np.int64)
    # values prod_x P[s,x]^Ex[a,x]
    V = np.ones((len(mons), len(S)), dtype=object)
    out = np.ones((len(mons), len(S)), dtype=np.int64)
    for x in range(P.shape[1]):
        col = P[:, x][None, :] ** Ex[:, x][:, None]   # N x |S|  (int64, small entries)
        out = (out * col) % p
    return out

def nmod(A, p):
    return flint.nmod_mat(A.shape[0], A.shape[1], [int(v) for v in A.flatten()], p)

def Me_modp(S, e, p, d):
    """nmod_mat whose columns span M_e(S) (rows = degree-e monomials)."""
    if len(S) == 0:
        return None
    lt = monos_lt(d, e)
    K = None
    if len(lt) == 0:
        nul = len(S)
        Kmat = flint.nmod_mat(len(S), nul, [1 if i == j else 0 for i in range(len(S)) for j in range(nul)], p)
    else:
        M = nmod(evalmat(S, lt, p), p)
        Kmat, nul = M.nullspace()
        if nul == 0:
            return None
        # keep first nul columns
        Kmat = flint.nmod_mat(len(S), nul, [int(Kmat[i, j]) for i in range(len(S)) for j in range(nul)], p)
    E = nmod(evalmat(S, monos(d, e), p), p)
    return E * Kmat

def colrank(Mm):
    return 0 if Mm is None else Mm.rank()

def hcat(mats, p, nrows):
    mats = [m for m in mats if m is not None]
    if not mats:
        return None
    cols = sum(m.ncols() for m in mats)
    ent = []
    for i in range(nrows):
        for m in mats:
            ent += [int(m[i, j]) for j in range(m.ncols())]
    return flint.nmod_mat(nrows, cols, ent, p)

def shift_dims(S, L, e, p, d):
    MS = Me_modp(S, e, p, d)
    hs = colrank(MS)
    nr = len(monos(d, e))
    # accumulate with periodic reduction to a basis (rank via transpose rref)
    acc = None; cur = 0
    for l in L:
        Ml = Me_modp(shiftset(S, l), e, p, d)
        if Ml is None:
            continue
        acc = hcat([acc, Ml], p, nr)
        # reduce: keep a column basis
        T = acc.transpose().rref()[0]
        r = acc.rank()
        acc = flint.nmod_mat(nr, r, [int(T[j, i]) for i in range(nr) for j in range(r)], p) if r > 0 else None
        cur = r
        if cur == hs:
            break
    return hs, cur
