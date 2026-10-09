"""Direct (non-equivariant) Hilbert function of a point set A in {sum k = s} mod p:
h_A(d) = rank of all monomials of degree <= d in k_2..k_M evaluated on A."""
import numpy as np
from modla import Echelon


def direct_profile(pts, p, batch=2000):
    A = np.array(pts, dtype=np.int64)
    n, M = A.shape
    X = A[:, 1:] % p
    E = Echelon(n, p)
    prof = []
    frontier = {(): np.ones(n, dtype=np.int64)}
    E.add(np.ones((1, n), dtype=np.int64))
    prof.append(E.rank)
    while E.rank < n:
        newf = {}
        for key, v in frontier.items():
            st = key[-1] if key else 0
            for j in range(st, M - 1):
                newf[key + (j,)] = (v * X[:, j]) % p
        frontier = newf
        vals = list(newf.values())
        for b in range(0, len(vals), batch):
            if E.rank == n:
                break
            E.add(np.array(vals[b:b + batch], dtype=np.int64))
        prof.append(E.rank)
    return prof
