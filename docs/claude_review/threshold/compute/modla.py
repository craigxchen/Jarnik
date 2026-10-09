"""Exact modular linear algebra over F_p, p < 2^31, with numpy.

mm(A, B, p): exact (A @ B) mod p.  Entries are split into 16-bit halves and multiplied with
float64 BLAS (Karatsuba: three products).  Every partial sum stays below 2^53, so each float
product is exact; the halves are recombined in int64.

Echelon: incremental reduced row echelon form.  add(B) inserts a batch of rows and returns the
indices of the rows of B that increased the rank (greedy, in row order), so the rank after a
batch is exact (mod p) and the set of returned rows is a basis of the span of all rows seen.
"""
import numpy as np

P1 = 2147483647          # 2^31 - 1
P2 = 2147483629          # 2^31 - 19

_CH_INNER = 1 << 17      # inner-dimension chunk: 2^17 * (2^17)^2 = 2^51 < 2^53
_OUT_LIMIT = 16_000_000   # output blocks of at most this many entries (memory)


def is_prime(n):
    if n < 2:
        return False
    small = (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37)
    for q in small:
        if n % q == 0:
            return n == q
    d, r = n - 1, 0
    while d % 2 == 0:
        d //= 2
        r += 1
    for a in small:
        x = pow(a, d, n)
        if x in (1, n - 1):
            continue
        for _ in range(r - 1):
            x = x * x % n
            if x == n - 1:
                break
        else:
            return False
    return True


assert is_prime(P1) and is_prime(P2)


def mm(A, B, p):
    """exact (A @ B) mod p for int64 arrays with entries in [0, p), p < 2^31."""
    assert p < (1 << 31)
    a_rows, k = A.shape
    k2, b_cols = B.shape
    assert k == k2
    out = np.zeros((a_rows, b_cols), dtype=np.int64)
    if k == 0 or a_rows == 0 or b_cols == 0:
        return out
    if a_rows * b_cols > _OUT_LIMIT and a_rows > 1:
        step = max(1, _OUT_LIMIT // b_cols)
        for r0 in range(0, a_rows, step):
            out[r0:r0 + step] = mm(A[r0:r0 + step], B, p)
        return out
    two32 = (1 << 32) % p
    for s in range(0, k, _CH_INNER):
        a = A[:, s:s + _CH_INNER]
        b = B[s:s + _CH_INNER]
        a1 = (a >> 16).astype(np.float64)
        a0 = (a & 0xFFFF).astype(np.float64)
        b1 = (b >> 16).astype(np.float64)
        b0 = (b & 0xFFFF).astype(np.float64)
        hh = (a1 @ b1).astype(np.int64)
        ll = (a0 @ b0).astype(np.int64)
        mid = ((a1 + a0) @ (b1 + b0)).astype(np.int64)
        cross = mid - hh - ll
        r = ((hh % p) * two32 + (cross % p) * 65536 + ll) % p
        out = (out + r) % p
    return out


def _split(A):
    a1 = (A >> 16).astype(np.float64)
    a0 = (A & 0xFFFF).astype(np.float64)
    return a1, a0, a1 + a0


def _combine(hh, ll, mid, p):
    hh = hh.astype(np.int64)
    ll = ll.astype(np.int64)
    cross = mid.astype(np.int64) - hh - ll
    return ((hh % p) * ((1 << 32) % p) + (cross % p) * 65536 + ll) % p


def _naive_chunk(C, p):
    """RREF of a small chunk C (rows) mod p.  Returns (R, pivots, idx): the nonzero
    normalised rows (mutually reduced), their pivot columns, and the indices of the chunk rows
    that produced them (greedy order)."""
    C = C.copy()
    c = C.shape[0]
    pivs, idx = [], []
    for i in range(c):
        row = C[i]
        nz = np.flatnonzero(row)
        if len(nz) == 0:
            continue
        j = int(nz[0])
        inv = pow(int(row[j]), p - 2, p)
        row = (row * inv) % p
        C[i] = row
        f = C[:, j].copy()
        f[i] = 0
        nzr = np.flatnonzero(f)
        if len(nzr):
            C[nzr] = (C[nzr] - (f[nzr, None] * row[None, :]) % p) % p
        pivs.append(j)
        idx.append(i)
    return C[idx], pivs, idx


def _rec_rref(B, p, room, base=48):
    """RREF of the rows of B (balanced recursion, early exit when rank == room).
    Returns (E, piv, idx): E rows in RREF on pivots piv, idx = indices of the rows of B that
    raised the rank in the greedy (row-order) sense."""
    m = B.shape[0]
    if m <= base:
        R, pivs, idx = _naive_chunk(B, p)
        if len(pivs) > room:
            R, pivs, idx = R[:room], pivs[:room], idx[:room]
            # rows beyond room cannot occur for a consistent room; kept for safety
        return R, np.array(pivs, dtype=np.int64), list(idx)
    h = m // 2
    E1, p1, i1 = _rec_rref(B[:h], p, room, base)
    if len(p1) >= room:
        return E1, p1, i1
    Bot = B[h:]
    if len(p1):
        X = Bot[:, p1]
        if X.any():
            Bot = (Bot - mm(X, E1, p)) % p
    if not Bot.any():
        return E1, p1, i1
    E2, p2, i2 = _rec_rref(Bot, p, room - len(p1), base)
    if len(p2) == 0:
        return E1, p1, i1
    if len(p1):
        X = E1[:, p2]
        if X.any():
            E1 = (E1 - mm(X, E2, p)) % p
        return np.vstack([E1, E2]), np.concatenate([p1, p2]), i1 + [h + i for i in i2]
    return E2, p2, [h + i for i in i2]


class Echelon:
    """row space of all rows added so far, kept as one matrix in RREF."""

    def __init__(self, ncols, p):
        assert p < (1 << 31)
        self.N = ncols
        self.p = p
        self.E = np.zeros((0, ncols), dtype=np.int64)
        self.piv = np.zeros(0, dtype=np.int64)

    @property
    def rank(self):
        return self.E.shape[0]

    def add(self, B):
        """insert rows of B; return indices of rows that raised the rank (greedy order)."""
        p = self.p
        B = np.asarray(B, dtype=np.int64) % p
        room = self.N - self.rank
        if B.shape[0] == 0 or room == 0:
            return []
        if self.rank:
            X = B[:, self.piv]
            if X.any():
                B = (B - mm(X, self.E, p)) % p
        if not B.any():
            return []
        EB, pB, idx = _rec_rref(B, p, room)
        if len(pB) == 0:
            return []
        if self.rank:
            X = self.E[:, pB]
            if X.any():
                self.E = (self.E - mm(X, EB, p)) % p
            self.E = np.vstack([self.E, EB])
            self.piv = np.concatenate([self.piv, pB])
        else:
            self.E, self.piv = EB, pB
        return idx

    def rref(self):
        return self.E, self.piv

    def kernel_basis(self):
        """basis K (rows) of {x : E x = 0}; shape (N - rank, N)."""
        p = self.p
        E, piv = self.E, self.piv
        N = self.N
        pivset = set(int(x) for x in piv)
        free = [j for j in range(N) if j not in pivset]
        K = np.zeros((len(free), N), dtype=np.int64)
        for a, f in enumerate(free):
            K[a, f] = 1
            K[a, piv] = (-E[:, f]) % p
        return K
