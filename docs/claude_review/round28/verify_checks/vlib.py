"""Round-28 verifier library (written from the definitions; imports nothing from round 27 or the primary lane).

Box convention: n points 0..n-1, base point n-1, d = n-1, k' in {-1,0,1}^d, kfull = (k', -sum k').
Profile = set of subsets S of {0..n-1}; box A = {k' : |kfull(S)| <= L for every S}.
"""
import itertools, random, sys, time
import numpy as np
import flint

Q1 = 2147483647   # 2^31 - 1
Q2 = 1073741827   # 2^30 + 3

# ---------------------------------------------------------------- profiles
def canon(n, S):
    S = tuple(sorted(S)); C = tuple(x for x in range(n) if x not in S)
    return min(S, C, key=lambda s: (len(s), s))

def syl(k):
    """Sylvester H_{2^k}: point x in F_2^k (binary of index), cut S_u = {x : <u,x> = 1}, u != 0."""
    n = 2 ** k
    pats = set()
    for u in range(1, n):
        pats.add(canon(n, [x for x in range(n) if bin(u & x).count('1') % 2 == 1]))
    return n, sorted(pats)

def leg(a, q):
    a %= q
    if a == 0: return 0
    return 1 if pow(a, (q - 1) // 2, q) == 1 else -1

def paleyK(q):
    """skew conference-type matrix K on (inf, 0..q-1): K = I + skew Paley core with border; K + K^T = 2I."""
    pts = ['inf'] + list(range(q)); n = q + 1
    K = np.zeros((n, n), dtype=int)
    for a, x in enumerate(pts):
        for b, y in enumerate(pts):
            if x == 'inf': K[a, b] = 1
            elif y == 'inf': K[a, b] = -1
            else: K[a, b] = (1 if x == y else 0) + leg(y - x, q)
    if not ((K + K.T) == 2 * np.eye(n, dtype=int)).all(): K = K.T
    assert ((K + K.T) == 2 * np.eye(n, dtype=int)).all()
    return K

def had(q):
    K = paleyK(q); n = q + 1
    H = K * K[0][None, :]
    assert (H[0] == 1).all() and ((H @ H.T) == n * np.eye(n, dtype=int)).all()
    return n, sorted(set(canon(n, [x for x in range(n) if H[r, x] == 1]) for r in range(1, n)))

def sdc(q):
    K = paleyK(q); n = q + 1
    D = np.diag(K[:, 0]); K = D @ K @ D
    pats = set()
    for a in range(1, n):
        S = [x for x in range(n) if K[x, a] == -1]
        pats.add(canon(n, S)); pats.add(canon(n, sorted(set(S) ^ {0, a})))
    return n, sorted(pats)

def restrict(n, pats, J):
    pos = {x: i for i, x in enumerate(J)}; m = len(J)
    out = set()
    for S in pats:
        T = [pos[x] for x in S if x in pos]
        if 0 < len(T) < m: out.add(canon(m, T))
    return m, sorted(out)

def rand_profile(n, seed):
    rng = random.Random(seed)
    r = rng.randint(n - 1, 2 * n)
    pats = set()
    while len(pats) < r:
        S = tuple(x for x in range(n) if rng.random() < 0.5)
        if 0 < len(S) < n: pats.add(canon(n, S))
    L = rng.choice([1, 2])
    return sorted(pats), L

# ---------------------------------------------------------------- boxes
def cube(d):
    return np.array(list(itertools.product([-1, 0, 1], repeat=d)), dtype=np.int8)

def box(n, pats, L):
    d = n - 1
    C = cube(d).astype(np.int16)
    full = np.concatenate([C, -C.sum(axis=1, keepdims=True)], axis=1)
    ok = np.ones(len(C), bool)
    for S in pats:
        ok &= np.abs(full[:, list(S)].sum(axis=1)) <= L
    return C[ok].astype(np.int8)

def spec_box(spec):
    """spec strings: syl:k:L, sdc:q:L, had:q:L, rsyl:k:J:L, rsdc:q:J:L, rhad:q:J:L (J comma list, 'inf'=0 idx), rand:n:seed"""
    f = spec.split(':')
    if f[0] == 'rand':
        n, seed = int(f[1]), int(f[2]); pats, L = rand_profile(n, seed)
        return n, pats, L, box(n, pats, L)
    if f[0] in ('syl', 'sdc', 'had'):
        n, pats = {'syl': syl, 'sdc': sdc, 'had': had}[f[0]](int(f[1])); L = int(f[2])
        return n, pats, L, box(n, pats, L)
    base = {'rsyl': syl, 'rsdc': sdc, 'rhad': had}[f[0]]
    N, P = base(int(f[1]))
    J = [int(x) for x in f[2].split(',')]
    n, pats = restrict(N, P, J); L = int(f[3])
    return n, pats, L, box(n, pats, L)

# ---------------------------------------------------------------- standard monomials
def std_monos(d, D):
    """exponent vectors with entries <= 2, degree <= D, grouped by degree (list per degree)."""
    by = [[] for _ in range(D + 1)]
    for a in itertools.product(range(3), repeat=d):
        s = sum(a)
        if s <= D: by[s].append(a)
    for b in by: b.sort(reverse=True)
    return by

def evalmat(S, monos):
    """|S| x len(monos) int8 matrix of k^alpha (k in {-1,0,1}: k^1 = k, k^2 = |k|)."""
    S = np.asarray(S, dtype=np.int8)
    absS = np.abs(S)
    out = np.empty((len(S), len(monos)), dtype=np.int8)
    for j, a in enumerate(monos):
        col = np.ones(len(S), dtype=np.int8)
        for x, ax in enumerate(a):
            if ax == 1: col = col * S[:, x]
            elif ax == 2: col = col * absS[:, x]
        out[:, j] = col
    return out

def to_nmod(M, p):
    vals = [0, 1, p - 1]
    flat = [vals[v] for v in M.ravel().tolist()]
    return flint.nmod_mat(M.shape[0], M.shape[1], flat, p)

def least_blocks(S, by, p, emax):
    """Return list over e=0..emax of numpy int64 arrays (H_S(e) x N_e) spanning M_e(S) (mod p), via one rref of the
    degree-ordered evaluation matrix: rows of rref with pivot in block e, restricted to block e."""
    monos = [a for e in range(emax + 1) for a in by[e]]
    starts = [0]
    for e in range(emax + 1): starts.append(starts[-1] + len(by[e]))
    E = evalmat(S, monos)
    R, r = to_nmod(E, p).rref()
    del E
    ncol = starts[-1]
    piv = []
    j = 0
    for i in range(r):
        while int(R[i, j]) == 0: j += 1
        piv.append(j); j += 1
    blocks = []
    for e in range(emax + 1):
        rows = [i for i in range(r) if starts[e] <= piv[i] < starts[e + 1]]
        B = np.zeros((len(rows), starts[e + 1] - starts[e]), dtype=np.int64)
        for a, i in enumerate(rows):
            for c in range(starts[e], starts[e + 1]):
                v = int(R[i, c])
                if v: B[a, c - starts[e]] = v
        blocks.append(B)
    return blocks

def shifted(A, l, index):
    """A_l = {k in A : k + l in A} as an array; index = dict tuple->i."""
    keep = []
    for i, k in enumerate(A):
        t = tuple((k.astype(int) + l).tolist())
        if t in index: keep.append(i)
    return A[keep]

def lpair_classes(d):
    out = []
    for x in range(d):
        v = [0] * d; v[x] = 1; out.append(tuple(v))
    for x in range(d):
        for y in range(x + 1, d):
            v = [0] * d; v[x] = 1; v[y] = -1; out.append(tuple(v))
    return out

def cube_classes(d):
    out = []
    for v in itertools.product([-1, 0, 1], repeat=d):
        if any(v):
            fnz = next(x for x in v if x)
            if fnz == 1: out.append(tuple(v))
    return out

def rank_rows(blocks, p):
    rows = [r for B in blocks for r in B.tolist()]
    if not rows: return 0
    return flint.nmod_mat(rows, p).rank()

def log(*a):
    print(time.strftime('%H:%M:%S'), *a, flush=True)
