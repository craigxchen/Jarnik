"""Round 28 compute lane: slab boxes, graded least spaces, multi-shift criterion (own pipeline).

Conventions (record: round 11/12, killer klib.build_box):
  n points, d = n-1, k' in {-1,0,1}^d, kfull = (k', -sum k'), slab p: |kfull(S_p)| <= L_p,
  S_p a canonical subset of {0..n-1} (the base point is n-1).
Least spaces: on the ternary cube k^3 = k, so Poly(S) is spanned by standard monomials (exponents <= 2) and the
degree-e moments of c in V(S,e) vanish on non-standard monomials.  One rref of the evaluation matrix
E (|S| x N_{<=E}, columns grouped by degree) gives, for every e, M_e(S) = span of the degree-e block of the rows
whose pivot lies in degree block e (row-echelon structure; see compute.md §2).
All linear algebra mod a prime p (flint nmod_mat).  M_e(A_{-l}) = M_e(A_l) since A_{-l} = A_l + l.
"""
import itertools, time, math, sys
import numpy as np
import flint

P1 = 2147483629
P2 = 2147483587

# ---------------------------------------------------------------- boxes
def canon(n, S):
    S = tuple(sorted(S))
    Sc = tuple(x for x in range(n) if x not in S)
    return min(S, Sc, key=lambda s: (len(s), s))

def fvec(n, S):
    """functional on k' (length d) giving kfull(S)"""
    d = n - 1
    f = [1 if x in S else 0 for x in range(d)]
    if (n - 1) in S:
        f = [v - 1 for v in f]
    return f

_cube = {}
def cube(d):
    if d not in _cube:
        _cube[d] = np.array(list(itertools.product((-1, 0, 1), repeat=d)), dtype=np.int8)
    return _cube[d]

def build_box(n, slabs):
    """slabs: list of (S, L). Returns int8 array of points (itertools.product order), and the functional data."""
    d = n - 1
    C = cube(d)
    if not slabs:
        return C.copy()
    F = np.array([fvec(n, S) for S, L in slabs], dtype=np.int16)
    Ls = np.array([L for S, L in slabs], dtype=np.int16)
    V = C.astype(np.int16) @ F.T
    mask = np.all(np.abs(V) <= Ls[None, :], axis=1)
    return C[mask].copy()

# ---------------------------------------------------------------- monomials
_mon = {}
def std_monos(d):
    """list (by degree e = 0..2d) of int8 arrays of exponent vectors in {0,1,2}^d"""
    if d not in _mon:
        allm = np.array(list(itertools.product((0, 1, 2), repeat=d)), dtype=np.int8)
        deg = allm.sum(axis=1)
        _mon[d] = [allm[deg == e] for e in range(2 * d + 1)]
    return _mon[d]

def Ncount(d, e):
    return len(std_monos(d)[e])

def evalmat(S, d, emax):
    """int8 matrix |S| x N_{<=emax} of standard monomials (degree-ordered blocks); also block boundaries."""
    mons = std_monos(d)
    Mx = np.concatenate(mons[:emax + 1], axis=0)
    bounds = np.cumsum([0] + [len(mons[e]) for e in range(emax + 1)])
    S = np.asarray(S, dtype=np.int8)
    V = np.ones((Mx.shape[0], S.shape[0]), dtype=np.int8)
    for x in range(d):
        T = np.stack([np.ones(S.shape[0], np.int8), S[:, x], S[:, x] * S[:, x]])  # 3 x |S|
        V *= T[Mx[:, x]]
    return V.T.copy(), bounds

def to_nmod(V, p):
    obj = np.array([0, 1, p - 1], dtype=object)[np.mod(V, 3)]
    return flint.nmod_mat(V.shape[0], V.shape[1], obj.ravel().tolist(), p)

def np_to_nmod(B, p):
    return flint.nmod_mat(B.shape[0], B.shape[1], [int(v) for v in B.ravel()], p)

def least_blocks(S, d, emax, p):
    """returns (H list for e=0..emax, blocks dict e -> int64 array H_e x N_e spanning M_e(S) mod p, rank)"""
    S = np.asarray(S, dtype=np.int8)
    if S.shape[0] == 0:
        return [0] * (emax + 1), {}, 0
    V, bounds = evalmat(S, d, emax)
    M = to_nmod(V, p)
    del V
    R, rk = M.rref()
    del M
    ncol = int(bounds[-1])
    piv = []
    j = 0
    for i in range(rk):
        while j < ncol and int(R[i, j]) == 0:
            j += 1
        piv.append(j)
        j += 1
    piv = np.array(piv, dtype=np.int64)
    H = []
    blocks = {}
    for e in range(emax + 1):
        a, b = int(bounds[e]), int(bounds[e + 1])
        rows = np.nonzero((piv >= a) & (piv < b))[0]
        H.append(len(rows))
        if len(rows):
            B = np.empty((len(rows), b - a), dtype=np.int64)
            for t, i in enumerate(rows):
                i = int(i)
                B[t, :] = [int(R[i, jj]) for jj in range(a, b)]
            blocks[e] = B
    return H, blocks, rk

def hilbert(S, d, p, emax=None):
    """graded Hilbert function H_S(e), e = 0..rho (full), via rref with growing column range"""
    S = np.asarray(S, dtype=np.int8)
    n = S.shape[0]
    if emax is None:
        e = 0
        while sum(Ncount(d, t) for t in range(e + 1)) < n and e < 2 * d:
            e += 1
    else:
        e = emax
    while True:
        H, blocks, rk = least_blocks(S, d, e, p)
        if rk == n or e >= 2 * d:
            break
        e += 1
    while H and H[-1] == 0:
        H.pop()
    return H, blocks

# ---------------------------------------------------------------- shifts
def Lpair(d):
    L = []
    for x in range(d):
        v = [0] * d; v[x] = 1; L.append(tuple(v))
    for x in range(d):
        for y in range(x + 1, d):
            v = [0] * d; v[x] = 1; v[y] = -1; L.append(tuple(v))
    return L   # one representative per class l ~ -l

def Lcube(d):
    out = []
    for v in itertools.product((-1, 0, 1), repeat=d):
        if any(v) and next(t for t in v if t) > 0:
            out.append(tuple(v))
    return out   # one representative per class l ~ -l

def shift_set(A, l):
    """A_l = {k in A : k + l in A} as int8 array"""
    A = np.asarray(A, dtype=np.int8)
    d = A.shape[1]
    key = lambda X: ((X.astype(np.int64) + 1) * (3 ** np.arange(d, dtype=np.int64))).sum(axis=1)
    kA = key(A)
    Ap = A.astype(np.int16) + np.array(l, dtype=np.int16)[None, :]
    ok = np.all(np.abs(Ap) <= 1, axis=1)
    kl = np.full(A.shape[0], -1, dtype=np.int64)
    kl[ok] = key(Ap[ok].astype(np.int8))
    sel = np.isin(kl, kA) & ok
    return A[sel]

class Acc:
    """per-degree accumulation of sum_l M_e(A_l) mod p, with lazy reduction"""
    def __init__(self, d, degs, target, p):
        self.d, self.p = d, p
        self.degs = list(degs)
        self.target = dict(target)
        self.basis = {e: None for e in self.degs}
        self.pend = {e: [] for e in self.degs}
        self.rank = {e: 0 for e in self.degs}
        self.full = {e: (self.target[e] == 0) for e in self.degs}
    def reduce(self, e):
        parts = ([self.basis[e]] if self.basis[e] is not None else []) + self.pend[e]
        self.pend[e] = []
        if not parts:
            return
        B = np.concatenate(parts, axis=0)
        M = np_to_nmod(B, self.p)
        R, rk = M.rref()
        Ne = B.shape[1]
        nb = np.empty((rk, Ne), dtype=np.int64)
        for i in range(rk):
            nb[i, :] = [int(R[i, j]) for j in range(Ne)]
        self.basis[e] = nb
        self.rank[e] = rk
        if rk >= self.target[e]:
            self.full[e] = True
    def add(self, blocks):
        for e in self.degs:
            if self.full[e] or e not in blocks:
                continue
            self.pend[e].append(blocks[e])
            npend = sum(b.shape[0] for b in self.pend[e])
            Ne = blocks[e].shape[1]
            if npend + self.rank[e] > max(2 * Ne, 2 * self.target[e]):
                self.reduce(e)
    def finish(self):
        for e in self.degs:
            if not self.full[e]:
                self.reduce(e)
    def allfull(self):
        return all(self.full[e] for e in self.degs)

def run_shifts(A, d, H, shifts, p, degs, acc=None, log=None, tag="", reduce_every=4):
    """add the shift classes (sorted by |A_l| desc) to the accumulator for degrees degs; returns acc"""
    if acc is None:
        acc = Acc(d, degs, {e: H[e] for e in degs}, p)
    emax = max(degs)
    sets = []
    for l in shifts:
        Al = shift_set(A, l)
        if Al.shape[0] > 0:
            sets.append((Al.shape[0], l, Al))
    sets.sort(key=lambda t: -t[0])
    t0 = time.time()
    for cnt, (sz, l, Al) in enumerate(sets):
        active = [e for e in degs if not acc.full[e] and sz > e]
        if not active:
            if acc.allfull():
                break
            continue
        hi = max(active)
        Hl, blocks, rk = least_blocks(Al, d, hi, p)
        acc.add(blocks)
        if (cnt + 1) % reduce_every == 0:
            for e in degs:
                if not acc.full[e] and acc.pend[e]:
                    acc.reduce(e)
        if log and ((cnt + 1) % 10 == 0 or cnt < 3):
            log("  %s shift %d/%d |A_l|=%d ranks %s (%.0fs)" % (tag, cnt + 1, len(sets), sz,
                {e: (acc.rank[e], acc.target[e]) for e in degs if not acc.full[e]}, time.time() - t0))
        if acc.allfull():
            break
    acc.finish()
    return acc

def reach_from(acc, degs, rho):
    for e in sorted(degs):
        if not acc.full[e]:
            return e
    return rho
