"""Exact (flint integer linear algebra): for every support R (|R| = s, rows 0..11, base row 0 included) the minimum,
over nonzero real z with supp z = R and sum z = 0, of t(z) = #{core columns c : <z, beta_c> != 0}.
Method: flats of the 12 restricted Hadamard columns h_c|_R containing the all-ones column; for a flat Z, generic z in
the orthogonal complement of span(h_c|_R : c in Z) is orthogonal exactly to the columns in Z; it has full support R
iff the complement lies in no coordinate hyperplane.  Every z is generic in the complement of the flat of columns it
is orthogonal to, so the minimum over z equals the minimum over flats.
usage: python3 hadflats2.py smin smax"""
import itertools, sys
import numpy as np
from flint import fmpz_mat
from fplib import pal_patterns
M_, pats = pal_patterns(11)
core = [P for P in pats if len(P) == 6]
H = np.array([[1]*12] + [[1 if x in P else -1 for x in range(12)] for P in core], dtype=np.int64)
assert (H @ H.T == 12 * np.eye(12, dtype=np.int64)).all()


def rk(rows):
    if not rows:
        return 0
    return fmpz_mat(len(rows), len(rows[0]), [x for r in rows for x in r]).rank()


def fullsupp_complement(rows, s):
    Mz = fmpz_mat(len(rows), s, [x for r in rows for x in r])
    X, nul = Mz.nullspace()
    if nul == 0:
        return False
    for i in range(s):
        if all(X[i, j] == 0 for j in range(nul)):
            return False
    return True


smin, smax = int(sys.argv[1]), int(sys.argv[2])
for s in range(smin, smax + 1):
    mn = 99
    ext = []
    hist = {}
    for R in itertools.combinations(range(12), s):
        V = [[int(x) for x in H[c, list(R)]] for c in range(12)]
        seen = set()
        stack = [[0]]
        localmax = -1
        while stack:
            Zl = stack.pop()
            r = rk([V[c] for c in Zl])
            if r >= s:
                continue
            Z = frozenset(c for c in range(12) if rk([V[x] for x in Zl] + [V[c]]) == r)
            if Z in seen:
                continue
            seen.add(Z)
            if len(Z) > localmax and fullsupp_complement([V[c] for c in sorted(Z)], s):
                localmax = len(Z)
            for c in range(12):
                if c not in Z:
                    stack.append(sorted(Z) + [c])
        if localmax >= 0:
            t = 12 - localmax
            hist[t] = hist.get(t, 0) + 1
            if t < mn:
                mn = t
                ext = [R]
            elif t == mn:
                ext.append(R)
    print('support', s, ': min t_core =', mn, ' #extremal supports', len(ext), ' histogram of min t over supports',
          sorted(hist.items()), flush=True)
    if s <= 8:
        print('   first extremal supports:', ext[:6], flush=True)
