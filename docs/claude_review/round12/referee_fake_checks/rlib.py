"""Referee's independent library (built from the doc definitions, not from fake_checks/).
Convention: n points 0..n-1, base point index `base`; k' = coordinates of the non-base points in {-1,0,1};
k_base = -sum k'.  A pattern S (subset of points) gives slab |k(S)| <= L (same for complement)."""
import itertools, numpy as np

def legendre(a, q):
    a %= q
    return 0 if a == 0 else (1 if pow(a, (q-1)//2, q) == 1 else -1)

def paley_I(q):
    # rows/cols indexed by inf (0) then F_q (1..q):  H(inf,.)=1, H(i,inf)=-1, H(i,j)=delta_ij+chi(j-i)
    M = q+1; H = np.zeros((M, M), dtype=int)
    for a in range(M):
        for b in range(M):
            if a == 0: H[a, b] = 1
            elif b == 0: H[a, b] = -1
            else:
                i, j = a-1, b-1
                H[a, b] = (1 if i == j else 0) + legendre(j-i, q)
    assert (H @ H.T == M*np.eye(M, dtype=int)).all()
    return H

def canon(n, S):
    S = frozenset(S); Sc = frozenset(range(n)) - S
    return min(S, Sc, key=lambda s: (len(s), sorted(s)))

def sdc(q):
    K = paley_I(q); M = q+1
    if not ((K + K.T) == 2*np.eye(M, dtype=int)).all(): K = K.T
    assert ((K + K.T) == 2*np.eye(M, dtype=int)).all()
    D = np.diag(K[:, 0]); K = D @ K @ D
    assert (K[:, 0] == 1).all() and ((K + K.T) == 2*np.eye(M, dtype=int)).all()
    pats = set()
    for a in range(1, M):
        S = {x for x in range(M) if K[x, a] == -1}
        Sp = S ^ {0, a}
        pats.add(canon(M, S)); pats.add(canon(M, Sp))
    assert len(pats) == 2*(M-1)
    return sorted(pats, key=lambda s: sorted(s))

def pal(q, astar_row_label=1, core_only=False):
    """construct.md Def: H_n = Paley-I with rows i in F_q negated; labels = columns 1..q;
    unflipped label columns + one flipped copy per row x (label a(x)=x for x!=inf, a(inf)=a*)."""
    H = paley_I(q); M = q+1
    Hn = H.copy(); Hn[1:, :] *= -1
    assert (Hn[:, 0] == 1).all()
    pats = set()
    for a in range(1, M):
        S = {x for x in range(M) if Hn[x, a] == 1}
        pats.add(canon(M, S))
    if not core_only:
        for x in range(M):
            lab = astar_row_label if x == 0 else x
            col = Hn[:, lab].copy(); col[x] *= -1
            S = {y for y in range(M) if col[y] == 1}
            if 0 < len(S) < M: pats.add(canon(M, S))
    return sorted(pats, key=lambda s: sorted(s))

def sylvester(k):
    n = 2**k
    pats = set()
    for a in range(1, n):
        S = {x for x in range(n) if bin(x & a).count('1') % 2 == 0}
        pats.add(canon(n, S))
    return sorted(pats, key=lambda s: sorted(s))

def box(n, pats, L, base):
    others = [x for x in range(n) if x != base]
    pm = np.zeros((len(pats), n), dtype=np.int64)
    for i, S in enumerate(pats):
        for x in S: pm[i, x] = 1
    K = np.array(list(itertools.product((-1, 0, 1), repeat=n-1)), dtype=np.int64)
    full = np.zeros((K.shape[0], n), dtype=np.int64)
    full[:, others] = K; full[:, base] = -K.sum(1)
    vals = full @ pm.T
    ok = (np.abs(vals) <= L).all(1)
    return K[ok], others

def reduced_monos(d, maxdeg):
    out = []
    for e in itertools.product((0, 1, 2), repeat=d):
        s = sum(e)
        if s <= maxdeg: out.append(e)
    out.sort(key=lambda e: (sum(e), e))
    return out

def eval_matrix(A, monos):
    A = np.asarray(A, dtype=np.int64)
    Mx = np.ones((len(monos), A.shape[0]), dtype=np.int64)
    for i, e in enumerate(monos):
        v = np.ones(A.shape[0], dtype=np.int64)
        for j, ej in enumerate(e):
            if ej: v = v * A[:, j]**ej
        Mx[i] = v
    return Mx
