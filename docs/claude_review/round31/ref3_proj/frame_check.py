# Frame bound sum_i <u,f_i>^2 <= 2|u|^2 when all off-diagonal G_ij < 0 (III.4/III.7):
# columns = all balanced cuts of m marks (plus random extra cuts), random positive weights.
import numpy as np, itertools
rng = np.random.default_rng(31)
worst = 0; tested = 0
for trial in range(400):
    m = int(rng.integers(4, 13))
    cols = []
    for S in itertools.combinations(range(m), m//2):
        c = -np.ones(m); c[list(S)] = 1; cols.append(c)
    extra = int(rng.integers(0, len(cols)//3+1))
    for _ in range(extra):
        c = np.where(rng.random(m) < 0.5, 1., -1.); cols.append(c)
    F = np.array(cols).T
    w = rng.random(F.shape[1])*float(rng.choice([0.3, 1.0, 3.0])) + 1; w /= w.sum()
    G = (F*w) @ F.T
    off = G - np.diag(np.diag(G))
    if (off[~np.eye(m, dtype=bool)] >= 0).any(): continue
    tested += 1; worst = max(worst, np.linalg.eigvalsh(G).max())
print("tested", tested, "worst top eigenvalue (must be <= 2):", worst)
