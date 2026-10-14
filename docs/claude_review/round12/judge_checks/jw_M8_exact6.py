"""Exact (over Z) dim of (J_w)_6 for M8: all degree-6 multiples of the width generators <k',v>^{h+1} (h+1 <= 6).
Compared with dim J(A)_6 = 637 (an exact value? J(A)_{<=6} cumulative: rank over Q of eval matrix computed exactly here too)."""
import itertools, math, sys
import numpy as np
from flint import fmpz_mat
sys.path.insert(0, '.')
from jw_test import widths, sylvester8, mono_tuple

def powlin(v, e, d):
    res = {tuple([0]*d): 1}
    for _ in range(e):
        new = {}
        for a, ca in res.items():
            for x in range(d):
                if v[x]:
                    k = list(a); k[x] += 1; k = tuple(k)
                    new[k] = new.get(k, 0) + ca*v[x]
        res = {k: c for k, c in new.items() if c}
    return res

A = sylvester8(2); d = 7; D = 6
gens = [(v, h) for v, h in widths(A, d, 2, D-1)]
m6 = [mono_tuple(c, d) for c in itertools.combinations_with_replacement(range(d), D)]
idx = {e: i for i, e in enumerate(m6)}
rows = []
for v, h in gens:
    g = powlin(v, h+1, d)
    for c in itertools.combinations_with_replacement(range(d), D-h-1):
        mm = mono_tuple(c, d)
        row = [0]*len(m6)
        for e, co in g.items():
            row[idx[tuple(p+q for p, q in zip(e, mm))]] += co
        rows.append(row)
M = fmpz_mat(len(rows), len(m6), [x for r in rows for x in r])
rkJw = M.rank()
# exact cumulative rank of evaluation matrix for degrees <=5 and <=6
Av = np.array(A, dtype=object)
def evalrank(D):
    cols = []
    for i in range(D+1):
        for c in itertools.combinations_with_replacement(range(d), i):
            col = [1]*len(A)
            for x in c: col = [col[t]*A[t][x] for t in range(len(A))]
            cols.append(col)
    E = fmpz_mat(len(cols), len(A), [x for c in cols for x in c])
    return E.rank()
r5, r6 = evalrank(5), evalrank(6)
H6 = r6 - r5
dimJA6 = len(m6) - H6
s = f"M8 degree 6, exact over Q: dim (J_w)_6 = {rkJw} (from {len(rows)} generator multiples); H_A(6) = {H6}; dim J(A)_6 = {dimJA6}; gap = {dimJA6 - rkJw}"
print(s); open('jw_M8_exact6_out.txt', 'w').write(s + "\n")
