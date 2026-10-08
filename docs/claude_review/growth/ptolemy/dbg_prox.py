import numpy as np, itertools
from prox_newton import run
rng = np.random.default_rng(1)
nr, sep, lsep, coef, allrho, rows = run(5, rng)
items = list(allrho.items())
for (T1, r1), (T2, r2) in itertools.combinations(items, 2):
    if abs(r1 - r2) < 1e-6:
        print(sorted(T1), sorted(T2), r1, r2)
