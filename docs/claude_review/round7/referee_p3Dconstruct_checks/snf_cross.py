# Cross-check of Z-surjectivity (index 1) of a half-size row subset via rank mod many primes
# and via sympy Smith normal form if available.
import random, sys
sys.path.insert(0, '.')
from paley_snf import paley, det_divisor_gcd, rank_mod
random.seed(7)
for q0 in [43, 47]:
    H = paley(q0); M = q0 + 1
    P = [[(1 + H[x][a])//2 for a in range(1, M)] for x in range(M)]
    order = list(range(M)); random.shuffle(order)
    R = []
    for x in order:
        if det_divisor_gcd([P[y] for y in R + [x]]) == 1: R.append(x)
    sub = [P[y] for y in R]
    ranks = {p: rank_mod(sub, p) for p in [2, 3, 5, 7, 11, 13, 23, 43, 47, 101]}
    try:
        from sympy import Matrix
        from sympy.matrices.normalforms import smith_normal_form
        S = smith_normal_form(Matrix(sub))
        diag = [S[i, i] for i in range(len(R))]
        sn = f"SNF diag product={eval('*'.join(str(abs(d)) for d in diag))}"
    except Exception as e:
        sn = f"sympy unavailable ({e.__class__.__name__})"
    print(f"q0={q0}: |R|={len(R)}, ranks mod p={ranks}, {sn}")
