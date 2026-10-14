"""Interpolation degree rho and Hilbert function of the parity box P_j = {y in [-j,j]^3 : y1+y2+y3 even}
(image of a Sylvester box under a Fano-line triple of columns when n is large; in M8 the image is P_2 minus 6 points).
Exact over Q (fmpz ranks of the cumulative evaluation matrices)."""
import itertools
from flint import fmpz_mat
def monos(i):
    return [(a, b, i - a - b) for a in range(i + 1) for b in range(i + 1 - a)]
out = []
for j in range(1, 5):
    Y = [y for y in itertools.product(range(-j, j + 1), repeat=3) if sum(y) % 2 == 0]
    cols = []; ranks = []; i = 0
    while True:
        cols += [[y[0]**e[0] * y[1]**e[1] * y[2]**e[2] for y in Y] for e in monos(i)]
        r = fmpz_mat(len(cols), len(Y), [x for c in cols for x in c]).rank()
        ranks.append(r)
        if r == len(Y): break
        i += 1
    H = [ranks[0]] + [ranks[t] - ranks[t - 1] for t in range(1, len(ranks))]
    s = f"j={j}: |P_j|={len(Y)} rho={len(ranks)-1} (full box: 6j={6*j}; 4j={4*j}) H={H}"
    print(s, flush=True); out.append(s)
open('parity_box_rho_out.txt', 'w').write("\n".join(out) + "\n")
