"""(a) direct Chern numbers vs ledger vs bound, and (IV.2) sides, on explicit
configurations with different contact patterns; (b) the family (c+u_i N, 1)."""
import sympy as sp
from front import t, full, report, contact_data, iv2_data, key
from engine import chern_side
from fractions import Fraction as Fr

out = []

# A1: graphs y = u_i - t^2/u_i  (columns (u_i^2 - t^2, u_i))
for us in ([1, 2, 3, 5, 7], [1, 2, 3, 5, 7, 11], [2, 3, 5, 7, 11, 13, 17, 19, 23]):
    cols = [(u * u - t ** 2, u) for u in us]
    cd, ch, ar = full(cols)
    report(f"A1 u_i - t^2/u_i m={len(us)}", ch, ar)

# A2: one 7-fold point at t=0 with nested groups 7 > 5 > 3 and a pair chain
rows = [(0, 0, 0, 1), (0, 0, 0, 2), (0, 0, 1, 0), (0, 1, 0, 0), (0, 1, 0, 5), (1, 0, 0, 0), (2, 0, 0, 0)]
cols = [(a * t + b * t ** 2 + c * t ** 3 + e * t ** 4, 1) for a, b, c, e in rows]
cd, ch, ar = full(cols)
report("A2 nested 7>5>3 at t=0 + chain", ch, ar)
# A2': a 9-fold point plus two 6-fold points
cols = [(t * (t - 1) * (t + 1) * (i + 1) + t ** 2 * (i % 3), 1) for i in range(9)]
cd, ch, ar = full(cols)
report("A2' 9 sections, 9-fold point at 0", ch, ar)

# A3: prime powers and tangencies
cols = [(u * (t - 1) ** 2 * (t + 2) ** 3, 1) for u in (1, 2, 3, 4)]
cols += [((t - 1) ** 2 * (t - 7), 1), ((t - 1) ** 2 * (t + 2) ** 3 + (t - 9) ** 3, 1)]
cd, ch, ar = full(cols)
report("A3 prime powers (t-1)^2 (t+2)^3, tangency order 3", ch, ar)

# A4: contact at infinity (a 3-group at the infinity point after normalisation, and a pair chain there)
fs = [t ** 3 + t, t ** 3 + t + 7, t ** 3 + 3 * t ** 2, 2 * t ** 3, 3 * t ** 3 + 1, -t ** 3 + t ** 2]
cols = [(f, 1) for f in fs]
cd, ch, ar = full(cols)
report("A4 contact at infinity", ch, ar)
# A4': non-polynomial columns, sections through several infinity points
cols = [(t ** 2 + 1, t), (t ** 2 + 2, 2 * t), (t ** 2 - 1, t + 3), (3 * t ** 2, t ** 2 + 1), (t ** 2 + t, t ** 2 - 5), (1, t ** 2 + 7)]
cd, ch, ar = full(cols)
report("A4' rational sections, mixed infinity", ch, ar)

# (b) columns (c + u_i N, 1), every row marked at every root (n_i = N)
for m in (6, 12, 30):
    Nn = (t - 1) * (t - 2) ** 2 * (t ** 2 + 1)
    c = t ** 3 + 5
    cols = [(c + (i + 1) * Nn, 1) for i in range(m)]
    cd = contact_data(cols)
    N = {}
    a = [dict() for _ in range(m)]
    for f, e in sp.Poly(Nn, t).factor_list()[1]:
        k = key(f)
        N[k] = (f.monic(), e)
        for i in range(m):
            a[i][k] = e
    ar = iv2_data(cd, N, a)
    ch = chern_side(m, cd['finite'], cd['vinf'], cd['degs'])
    report(f"B (c+u_iN,1) m={m}", ch, ar)
    print("   predicted slack 18L =", 18 * ar['L'], " lhs-leading m^2L/4 =", Fr(m * m * ar['L'], 4),
          " rhs-leading m^2 delta =", m * m * ar['delta'])
