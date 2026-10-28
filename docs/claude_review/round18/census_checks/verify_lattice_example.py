"""Independent exact check of one lattice-census output at k = 8 (Data): rebuild the best vector for
one random prime-block vector, verify n_T | X_i - X_j for all i, j in T exactly, verify the
datum identity det(P_i,P_j) = t_ij * prod_{T contains i,j} n_T with integer t_ij, and print H, w, eta."""
import random
from itertools import product
from math import log
from sympy import primerange
from lattice_census import visible, lattice_basis
k, A, seed = 8, 1500, 7
m = k - 1
rng = random.Random(seed)
Ts = visible(m)
chosen = rng.sample(list(primerange(A, 2 * A)), len(Ts))
blocks = dict(zip(Ts, chosen))
basis = lattice_basis(m, blocks)
best = None
for c in product(range(-2, 3), repeat=m - 1):
    if not any(c):
        continue
    X = [0] + [sum(c[l] * basis[l][r] for l in range(m - 1)) for r in range(m - 1)]
    if len(set(X)) < m:
        continue
    H = 0.0
    ok = True
    for i in range(m):
        for j in range(i + 1, m):
            G = 1
            for T, p in blocks.items():
                if (i + 1) in T and (j + 1) in T:
                    G *= p
            if (X[i] - X[j]) % G:
                ok = False
            H += log(abs(X[i] - X[j]) // G)
    assert ok
    if best is None or H < best[0]:
        best = (H, X)
H, X = best
w = log(A)
eta = max(abs(log(p) - w) for p in chosen) / w
print("k=8 A=%d: verified divisibility for best vector; H_Sigma = %.2f, w = %.3f, eta = %.3f, H/w = %.1f" % (A, H, w, eta, H / w))
print("max log10|X_i| = %.1f" % (max(log(abs(x)) for x in X if x) / log(10)))
