"""Killer's independent library (written from scratch; imports nothing from interp_checks/).

Conventions: a finite set S of points of Z^d (tuples).  For degree e,
  V(S,e)  = measures c on S with all moments of order < e vanishing,
  M_e(S)  = { mu_e(c) : c in V(S,e) }   (degree-e moment vectors; = least space in degree e, as coefficient
            vectors w.r.t. the monomial basis of degree-e moments),
  H_S(e)  = dim M_e(S)   (graded Hilbert function of the point set).
(SHIFT_e) w.r.t. shifts L: M_e(S) = sum_{l in L} M_e(S cap (S - l)).
Linear algebra: exact over Q (fmpz nullspace + fmpq rank) or mod a prime p (nmod_mat).
"""
import itertools, math
import flint

P1 = 2147483647   # 2^31 - 1
P2 = 1000000007

def monos(d, e):
    out = []
    for comb in itertools.combinations_with_replacement(range(d), e):
        a = [0] * d
        for i in comb:
            a[i] += 1
        out.append(tuple(a))
    return out

def monos_lt(d, e):
    out = []
    for t in range(e):
        out += monos(d, t)
    return out

def mval(k, a):
    v = 1
    for ki, ai in zip(k, a):
        if ai:
            v *= ki ** ai
    return v

def shiftset(S, l):
    """S_l = S cap (S - l) = {k in S : k + l in S}"""
    Ss = set(S)
    return [k for k in S if tuple(ki + li for ki, li in zip(k, l)) in Ss]

# ---------------- mod p ----------------
def _nullspace_modp(rows, ncols, p):
    """rows: list of lists (r x ncols). Returns list of basis vectors of the right kernel (each length ncols)."""
    if not rows:
        return [[1 if i == j else 0 for i in range(ncols)] for j in range(ncols)]
    M = flint.nmod_mat(rows, p)
    X, nul = M.nullspace()
    return [[int(X[i, j]) for i in range(ncols)] for j in range(nul)]

def Me_basis_modp(S, e, p=P1):
    """basis (as list of vectors in F_p^{#monos(d,e)}) spanning M_e(S) mod p."""
    if not S:
        return []
    d = len(S[0])
    lt = monos_lt(d, e)
    rows = [[mval(k, a) % p for k in S] for a in lt]
    K = _nullspace_modp(rows, len(S), p)
    if not K:
        return []
    me = monos(d, e)
    E = [[mval(k, a) % p for k in S] for a in me]
    # columns: E * K_j
    vecs = []
    for kv in K:
        vecs.append([sum(E[i][t] * kv[t] for t in range(len(S)) if kv[t]) % p for i in range(len(me))])
    return vecs

def rank_modp(vecs, p=P1):
    if not vecs:
        return 0
    return flint.nmod_mat(vecs, p).rank()

def H_modp(S, e, p=P1):
    return rank_modp(Me_basis_modp(S, e, p), p)

def shift_dims_modp(S, L, e, p=P1):
    """returns (H_S(e), dim sum_l M_e(S_l)) mod p"""
    hs = H_modp(S, e, p)
    allv = []
    for l in L:
        Sl = shiftset(S, l)
        allv += Me_basis_modp(Sl, e, p)
    return hs, rank_modp(allv, p)

# ---------------- exact over Q ----------------
def _nullspace_Q(rows, ncols):
    if not rows:
        return [[1 if i == j else 0 for i in range(ncols)] for j in range(ncols)]
    M = flint.fmpz_mat(rows)
    X, nul = M.nullspace()
    return [[int(X[i, j]) for i in range(ncols)] for j in range(nul)]

def Me_basis_Q(S, e):
    if not S:
        return []
    d = len(S[0])
    lt = monos_lt(d, e)
    rows = [[mval(k, a) for k in S] for a in lt]
    K = _nullspace_Q(rows, len(S))
    if not K:
        return []
    me = monos(d, e)
    E = [[mval(k, a) for k in S] for a in me]
    return [[sum(E[i][t] * kv[t] for t in range(len(S)) if kv[t]) for i in range(len(me))] for kv in K]

def rank_Q(vecs):
    if not vecs:
        return 0
    return flint.fmpz_mat(vecs).rank()

def H_Q(S, e):
    return rank_Q(Me_basis_Q(S, e))

def shift_dims_Q(S, L, e):
    hs = H_Q(S, e)
    allv = []
    for l in L:
        allv += Me_basis_Q(shiftset(S, l), e)
    return hs, rank_Q(allv)

def hilbert(S, emax, exact=False, p=P1):
    return [H_Q(S, e) if exact else H_modp(S, e, p) for e in range(emax + 1)]

# ---------------- slab boxes in the round-11/12 convention ----------------
def kfull(kp):
    return tuple(kp) + (-sum(kp),)

def kS(kp, S):
    kf = kfull(kp)
    return sum(kf[x] for x in S)

def build_box(n, slabs, T=1):
    d = n - 1
    return [kp for kp in itertools.product(range(-T, T + 1), repeat=d)
            if all(abs(kS(kp, S)) <= L for S, L in slabs)]

def Lpair(d):
    L = []
    for x in range(d):
        e = [0] * d; e[x] = 1; L.append(tuple(e)); L.append(tuple(-v for v in e))
    for x in range(d):
        for y in range(d):
            if x != y:
                e = [0] * d; e[x] = 1; e[y] = -1; L.append(tuple(e))
    return L

def Lcoord(d):
    L = []
    for x in range(d):
        e = [0] * d; e[x] = 1; L.append(tuple(e)); L.append(tuple(-v for v in e))
    return L
