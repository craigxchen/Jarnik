import numpy as np, itertools
from balanced_newton import solve
rng = np.random.default_rng(5)
for tr in range(6):
    nr, sep, lsep, coef, allrho, rows = solve(5, rng)
    coll = []
    items = list(allrho.items())
    for (T1, r1), (T2, r2) in itertools.combinations(items, 2):
        if abs(r1 - r2) < 1e-6:
            coll.append((''.join(map(str, sorted(T1))), ''.join(map(str, sorted(T2)))))
    print(f"res {nr:.1e} lsep {lsep:.2e} collisions: {coll}")
