"""Full polynomial ring Q[y_1..y_d] graded pieces up to maxdeg+1 (no truncation), same interface as loclib.Ring."""
import itertools, numpy as np
from loclib import PR
class FullRing:
    def __init__(self, d, maxdeg):
        self.d = d; self.maxdeg = maxdeg; B = maxdeg + 2; self.B = B
        self.pw = B ** np.arange(d, dtype=np.int64)
        self.mons = []; self.idx = []
        for t in range(maxdeg + 2):
            arr = np.array([tuple(sum(1 for z in c if z == q) for q in range(d)) for c in itertools.combinations_with_replacement(range(d), t)], dtype=np.int64).reshape(-1, d)
            self.mons.append(arr)
            self.idx.append({int(c): i for i, c in enumerate(arr @ self.pw)})
        self.N = [len(m) for m in self.mons]
        self.mul = []
        for t in range(maxdeg + 1):
            tab = np.zeros((d, self.N[t]), dtype=np.int64)
            codes = self.mons[t] @ self.pw
            for l in range(d):
                tab[l] = [self.idx[t + 1][int(c + self.pw[l])] for c in codes]
            self.mul.append(tab)
    def one(self, mod=PR): return np.ones(1, dtype=np.int64)
    def mul_var(self, vec, t, l, mod=PR):
        out = np.zeros(self.N[t + 1], dtype=np.int64); out[self.mul[t][l]] = vec % mod; return out
    def mul_lin(self, vec, t, v, mod=PR):
        out = np.zeros(self.N[t + 1], dtype=np.int64)
        for l in range(self.d):
            if v[l] % mod:
                out[self.mul[t][l]] = (out[self.mul[t][l]] + vec * (v[l] % mod)) % mod
        return out
