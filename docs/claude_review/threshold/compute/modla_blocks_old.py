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


class Echelon:
    """row space of all rows added so far, stored as blocks.  Block j is in RREF on its own
    pivots and vanishes on the pivots of blocks < j.  Reducing a row through the blocks in
    order therefore clears all pivots.  Blocks are merged (fully reduced) when there are many."""

    def __init__(self, ncols, p, chunk=32, maxblocks=12):
        assert p < (1 << 31)
        self.N = ncols
        self.p = p
        self.blocks = []          # list of (E, piv, (f1, f0, fs))
        self.chunk = chunk
        self.maxblocks = maxblocks

    @property
    def rank(self):
        return sum(b[0].shape[0] for b in self.blocks)

    @property
    def piv(self):
        if not self.blocks:
            return np.zeros(0, dtype=np.int64)
        return np.concatenate([b[1] for b in self.blocks])

    def _reduce(self, B, start=0):
        p = self.p
        for E, pv, (f1, f0, fs) in self.blocks[start:]:
            X = B[:, pv]
            if not X.any():
                continue
            x1, x0, xs = _split(X)
            B = (B - _combine(x1 @ f1, x0 @ f0, xs @ fs, p)) % p
        return B

    def _merge_two(self, i):
        """merge blocks i and i+1 (adjacent) into one block, fully reduced."""
        p = self.p
        Ei, pi, _ = self.blocks[i]
        Ej, pj, sj = self.blocks[i + 1]
        X = Ei[:, pj]
        if X.any():
            x1, x0, xs = _split(X)
            Ei = (Ei - _combine(x1 @ sj[0], x0 @ sj[1], xs @ sj[2], p)) % p
        E = np.vstack([Ei, Ej])
        pv = np.concatenate([pi, pj])
        self.blocks[i:i + 2] = [(E, pv, _split(E))]

    def _geometric(self, start=0):
        while len(self.blocks) - start >= 2:
            a = self.blocks[-2][0].shape[0]
            b = self.blocks[-1][0].shape[0]
            if 2 * b >= a or len(self.blocks) - start > self.maxblocks:
                self._merge_two(len(self.blocks) - 2)
            else:
                break

    def _merge(self):
        while len(self.blocks) > 1:
            self._merge_two(len(self.blocks) - 2)

    def add(self, B):
        """insert rows of B; return indices of rows that raised the rank (greedy order)."""
        p = self.p
        B = np.asarray(B, dtype=np.int64) % p
        if B.shape[0] == 0 or self.rank == self.N:
            return []
        B = self._reduce(B)
        start = len(self.blocks)
        indep = []
        for cs in range(0, B.shape[0], self.chunk):
            if self.rank == self.N:
                break
            C = B[cs:cs + self.chunk]
            if not C.any():
                continue
            C = self._reduce(C, start)
            R, pivs, idx = _naive_chunk(C, p)
            if len(pivs):
                pv = np.array(pivs, dtype=np.int64)
                self.blocks.append((R, pv, _split(R)))
                self._geometric(start)
                indep.extend(cs + i for i in idx)
        self._geometric(0)
        return indep

    def rref(self):
        """(E, piv) fully reduced."""
        if not self.blocks:
            return np.zeros((0, self.N), dtype=np.int64), np.zeros(0, dtype=np.int64)
        if len(self.blocks) > 1:
            self._merge()
        return self.blocks[0][0], self.blocks[0][1]

    def kernel_basis(self):
        """basis K (rows) of {x : E x = 0}; shape (N - rank, N)."""
        p = self.p
        E, piv = self.rref()
        N = self.N
        pivset = set(int(x) for x in piv)
        free = [j for j in range(N) if j not in pivset]
        K = np.zeros((len(free), N), dtype=np.int64)
        for a, f in enumerate(free):
            K[a, f] = 1
            K[a, piv] = (-E[:, f]) % p
        return K
