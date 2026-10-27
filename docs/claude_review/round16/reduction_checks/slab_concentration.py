"""Exact: concentration of the slab polynomial Q_S(k) = prod_{|t|<=j} (k(S) - t) on the uniform ternary cube,
at regime parameters (T = 1, w = sqrt(8 d ln 4r), j = floor(w/2), m' = ceil(4w)), |S| = s.
Pr(Q_S = 0) = Pr(|k(S)| <= j) computed exactly (trinomial coefficients, big integers).
Shows: a polynomial of degree 2j+1 < m' vanishes on all but a tiny fraction of the cube, so no anti-concentration
bound for the cube measure that depends only on (degree, variance, influences) can be < 1 - Pr(Q_S != 0) at degree 2j+1."""
import math, sys
from fractions import Fraction

def trinom_row(s):
    # coefficients of (1 + t + t^2)^s, index i <-> sum = i - s
    row = [1]
    for _ in range(s):
        new = [0] * (len(row) + 2)
        for i, c in enumerate(row):
            new[i] += c; new[i + 1] += c; new[i + 2] += c
        row = new
    return row

out = []
def log(x):
    print(x); out.append(x); sys.stdout.flush()

for n, rlab in [(2901, 'n'), (5000, 'n'), (6000, 'n^2')]:
    d = n - 1
    r = n if rlab == 'n' else n * n
    w = math.sqrt(8 * d * math.log(4 * r))
    j = int(w // 2)
    mp = math.ceil(4 * w)
    for s in (d // 2, d):
        row = trinom_row(s)
        tot = 3 ** s
        inside = sum(row[s - j: s + j + 1])
        p_out = Fraction(tot - inside, tot)
        lg = math.log10(tot - inside) - s * math.log10(3)
        log(f"n={n} r={rlab} d={d} w={w:.2f} j={j} m'={mp} deg Q=2j+1={2*j+1}  s={s}:  Pr(Q_S != 0) = 10^{lg:.3f}"
            f"  (Hoeffding bound 2exp(-(j+1)^2/(2s)) = 10^{math.log10(2)-(j+1)**2/(2*s)/math.log(10):.3f})")
open(__file__.replace('.py', '_out.txt'), 'w').write('\n'.join(out) + '\n')
