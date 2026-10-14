# Incremental mod-p spanning certificate per isotypic component (cyclic sigma x {+-1}).
# Certificate: rank of the projected monomial columns == number of valid orbit rows  => U_chi = F^A_chi.
import numpy as np, itertools, time
from math import comb
from cyc_rank import is_prime, prim_root_of_unity
def find_prime(q, below=700000):
    p = below - (below % (2*q)) + 1
    while not is_prime(p): p -= 2*q
    return p
class Echelon:
    def __init__(self, R, p):
        self.R = R; self.p = p; self.B = np.zeros((0, R)); self.piv = []
    def add(self, X):            # X: R x b float64 (entries in [0,p))
        p = self.p
        if self.piv:
            P = np.array(self.piv)
            X = (X - (self.B.T @ X[P, :]) % p) % p
        Y = X.T.copy()           # b x R
        newrows = []; newpiv = []
        for i in range(Y.shape[0]):
            y = Y[i] % p
            for r_, pc in zip(newrows, newpiv):   # reduce against new rows so far
                if y[pc]: y = (y - y[pc] * r_) % p
            nzs = np.nonzero(y)[0]
            if len(nzs) == 0: continue
            pc = nzs[0]; inv = pow(int(y[pc]), p - 2, p)
            y = (y * inv) % p
            for k in range(len(newrows)):
                if newrows[k][pc]: newrows[k] = (newrows[k] - newrows[k][pc] * y) % p
            newrows.append(y); newpiv.append(pc)
        if newrows:
            N = np.array(newrows)
            if self.piv:
                self.B = (self.B - (self.B[:, newpiv] @ N) % p) % p
            self.B = np.vstack([self.B, N]); self.piv += newpiv
        return len(newrows)
