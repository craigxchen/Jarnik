import itertools, numpy as np, models
from cyc_run import fam_orbit_reps
for L in (2,):
    d, cols, Ls = models.paley12(L); A = models.box(d, cols, Ls)
    Y = A @ np.array(cols).T
    print("Paley L=2 |A|", len(A), " columns basis? det", round(np.linalg.det(np.array(cols, dtype=float))))
    for b in range(2, 12):
        stats = {}
        for S in fam_orbit_reps(cols, b, 11):
            P = set(map(tuple, Y[:, list(S)].tolist()))
            stats[len(P)] = stats.get(len(P), 0) + 1
        print(f"  b={b}: full box size {5**b}; image sizes -> #orbits {stats}", flush=True)
