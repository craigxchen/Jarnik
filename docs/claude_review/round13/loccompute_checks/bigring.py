"""Ring for large d (d up to ~16): reduced monomials (exponents 0..2) by degree, numpy code->index map.
Same interface as loclib.Ring (mons, N, mul, pw, idx replaced by code2pos array, one, mul_lin)."""
import numpy as np
from loclib import PR

def mons_by_degree(d, maxdeg):
    # dynamic programming over variables
    out = [[np.zeros((1, 0), dtype=np.int8)] + [np.zeros((0, 0), dtype=np.int8)] * maxdeg]
    cur = {0: np.zeros((1, 0), dtype=np.int8)}
    for x in range(d):
        new = {}
        for t, arr in cur.items():
            for e in (0, 1, 2):
                if t + e <= maxdeg:
                    a2 = np.hstack([arr, np.full((len(arr), 1), e, dtype=np.int8)])
                    new.setdefault(t + e, []).append(a2)
        cur = {t: np.vstack(v) for t, v in new.items()}
    return [cur.get(t, np.zeros((0, d), dtype=np.int8)) for t in range(maxdeg + 1)]

class BigRing:
    def __init__(self, d, maxdeg):
        self.d = d; self.maxdeg = maxdeg
        self.pw = 3 ** np.arange(d, dtype=np.int64)
        self.mons = [m.astype(np.int64) for m in mons_by_degree(d, maxdeg + 1)]
        self.N = [len(m) for m in self.mons]
        self.pos = np.full(3 ** d, -1, dtype=np.int32)
        self.codes = []
        for t in range(maxdeg + 2):
            c = self.mons[t] @ self.pw
            self.codes.append(c)
            self.pos[c] = np.arange(len(c), dtype=np.int32)
        self.mul = []
        for t in range(maxdeg + 1):
            tab = np.full((d, self.N[t]), -1, dtype=np.int64)
            for l in range(d):
                ok = self.mons[t][:, l] < 2
                tab[l, ok] = self.pos[self.codes[t][ok] + self.pw[l]]
            self.mul.append(tab)
        self.idx = None
    def index_of(self, t, e):
        return int(self.pos[int(np.asarray(e) @ self.pw)])
    def one(self, mod=PR):
        return np.ones(1, dtype=np.int64)
    def mul_lin(self, vec, t, v, mod=PR):
        out = np.zeros(self.N[t + 1], dtype=np.int64)
        for l in range(self.d):
            if v[l] % mod:
                tab = self.mul[t][l]; ok = tab >= 0
                out[tab[ok]] += (vec[ok] * (v[l] % mod)) % mod   # tab injective on ok
                out %= mod
        return out
