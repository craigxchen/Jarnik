import itertools, numpy as np, models
from sympy import Matrix
from sympy.matrices.normalforms import smith_normal_form
from cyc_run import fam_orbit_reps
d, cols, Ls = models.paley12(2); A = models.box(d, cols, Ls)
Y = A @ np.array(cols).T
print("columns (i+D):", [ [x for x in range(11) if c[x]] for c in cols])
for S in fam_orbit_reps(cols, 5, 11):
    P = np.array(sorted(set(map(tuple, Y[:, list(S)].tolist()))))
    if len(P) < 2000:
        diffs = Matrix((P[1:] - P[0]).tolist()[:200])
        snf = smith_normal_form(diffs.T * 1)
        print(S, len(P), "SNF diag of lattice spanned by P-P:", [snf[i,i] for i in range(5)])
        # find relation mod 3
        for c in itertools.product(range(3), repeat=5):
            if any(c) and all((np.dot(p, c) - np.dot(P[0], c)) % 3 == 0 for p in P):
                print("    relation mod 3: coeffs", c, " value", np.dot(P[0], c) % 3); break
# sizes of x-level: the full-lattice index
M = Matrix(cols)
print("SNF of column matrix:", [smith_normal_form(M)[i,i] for i in range(11)])
