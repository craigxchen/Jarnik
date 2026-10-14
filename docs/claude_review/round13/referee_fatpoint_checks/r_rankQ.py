"""Referee: exact (fmpz) rank of the evaluation matrix of grid-reduced monomials of degree <= e on Y_3(2)."""
import itertools, numpy as np, pickle, sys, time
from flint import fmpz_mat, nmod_mat
Y, D, T = pickle.load(open(sys.argv[1], 'rb'))
e = int(sys.argv[2])
mons = [a for a in itertools.product(range(5), repeat=7) if sum(a) <= e]
Yo = Y.astype(object)
rows = []
for a in mons:
    col = np.ones(len(Y), dtype=object)
    for x, ax in enumerate(a):
        if ax: col = col * Yo[:, x] ** ax
    rows.append([int(v) for v in col])
print('rows', len(rows), 'cols', len(Y), flush=True)
for p in (2147475331, 1000000007):
    t0 = time.time(); r = nmod_mat(rows, p).rank(); print('rank mod', p, '=', r, '%.0fs' % (time.time() - t0), flush=True)
if len(sys.argv) > 3:
    t0 = time.time(); r = fmpz_mat(rows).rank(); print('rank over Z =', r, '%.0fs' % (time.time() - t0), flush=True)
