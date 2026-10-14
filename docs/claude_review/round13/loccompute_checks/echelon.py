import numpy as np
from loclib import PR
class Echelon:
    """incremental row echelon basis mod PR (numpy int64). add(v) -> True if v independent (then stored)."""
    def __init__(self, n):
        self.n = n; self.rows = np.zeros((0, n), dtype=np.int64); self.piv = []
    def reduce(self, v):
        v = np.asarray(v, dtype=np.int64) % PR
        for r, p in zip(self.rows, self.piv):
            c = v[p]
            if c: v = (v - c * r) % PR
        return v
    def add(self, v):
        v = self.reduce(v)
        nz = np.nonzero(v)[0]
        if len(nz) == 0: return False
        p = nz[0]; inv = pow(int(v[p]), PR - 2, PR); v = (v * inv) % PR
        self.rows = np.vstack([self.rows, v[None, :]]); self.piv.append(p)
        return True
    @property
    def rank(self): return len(self.piv)
