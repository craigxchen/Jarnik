"""Library: least-space / leading-form computations mod a prime, and the multi-shift criterion.
M_e(S) = { (sum_k c_k k^alpha)_{|alpha|=e} : c in V(S,e) },  V(S,e) = {c : sum_k c_k k^beta = 0, |beta|<e}.
(SHIFT_e): sum_{l in L} M_e(S cap (S-l)) == M_e(S).
All ranks mod P (a large prime). Rank mod P <= rank over Q; see notes for which direction is rigorous."""
import itertools, sys
from flint import nmod_mat

P = 2147483629  # prime < 2^31

def monos_deg(d, e):
    """exponent tuples of total degree exactly e in d variables"""
    if d == 1:
        return [(e,)]
    out = []
    for a in range(e, -1, -1):
        for rest in monos_deg(d - 1, e - a):
            out.append((a,) + rest)
    return out

def monos_upto(d, e):  # degree < e
    out = []
    for t in range(e):
        out += monos_deg(d, t)
    return out

def powtab(S):
    return S

def evalmat(S, monos):
    """numpy evaluation; exact for integer points (values reduced mod P)."""
    import numpy as np
    big = any(abs(v) > 1 for k in S for v in k)
    K = np.array(S, dtype=object) if big else np.array(S, dtype=np.int64)
    out = np.ones((len(S), len(monos)), dtype=K.dtype)
    for j, a in enumerate(monos):
        col = np.ones(len(S), dtype=K.dtype)
        for x, ax in enumerate(a):
            if ax:
                col = col * (K[:, x] ** ax)
        out[:, j] = col
    return [[int(v) % P for v in row] for row in out.tolist()]

def mat(rows, nr, nc):
    if nr == 0 or nc == 0:
        return nmod_mat(max(nr,0), max(nc,0), P)
    return nmod_mat(nr, nc, [x for r in rows for x in r], P)

def M_space(S, e, d, cache=None):
    """returns nmod_mat whose rows span M_e(S) (mod P), possibly with dependent rows"""
    S = list(S)
    if len(S) == 0:
        return None
    Ne = monos_deg(d, e)
    Ee = evalmat(S, Ne)
    if e == 0:
        return mat([[1]], 1, 1)
    Nlow = monos_upto(d, e)
    El = evalmat(S, Nlow)
    # V(S,e) = right nullspace of El^T  (El^T : |Nlow| x |S|)
    ElT = mat([[El[i][j] for i in range(len(S))] for j in range(len(Nlow))], len(Nlow), len(S))
    X, nul = ElT.nullspace()
    if nul == 0:
        return None
    # basis vectors are the first nul columns of X (|S| x |S|)
    Xc = nmod_mat(nul, len(S), [int(X[i, j]) for j in range(nul) for i in range(len(S))], P)
    EeM = mat(Ee, len(S), len(Ne))
    return Xc * EeM   # nul x |Ne|

def rank_of(m):
    return 0 if m is None else m.rank()

def stack(mats, ncols):
    rows = []
    for m in mats:
        if m is None:
            continue
        for i in range(m.nrows()):
            rows.append([int(m[i, j]) for j in range(m.ncols())])
    if not rows:
        return None
    return nmod_mat(len(rows), ncols, [x for r in rows for x in r], P)

def shift_test(A, e, d, L, verbose=False):
    """returns (dim M_e(A), dim sum_l M_e(A_l))"""
    Aset = set(A)
    MA = M_space(A, e, d)
    rA = rank_of(MA)
    Ne = len(monos_deg(d, e))
    acc = None
    racc = 0
    for l in L:
        Al = [k for k in A if tuple(ki + li for ki, li in zip(k, l)) in Aset]
        if len(Al) == 0:
            continue
        Ml = M_space(Al, e, d)
        if Ml is None:
            continue
        acc = stack([acc, Ml], Ne) if acc is not None else Ml
        # reduce acc to a basis to keep it small
        R = acc.rref()[0]
        rr = acc.rank()
        acc = nmod_mat(rr, Ne, [int(R[i, j]) for i in range(rr) for j in range(Ne)], P) if rr > 0 else None
        racc = rr
        if racc == rA:
            break
    return rA, racc

def HF(S, e, d):
    """rank of Poly_{<=e} on S"""
    mon = monos_upto(d, e + 1)
    E = evalmat(S, mon)
    return mat(E, len(S), len(mon)).rank()
