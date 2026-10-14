"""Exact list of all flats giving exceptional generators (t_core <= 5) with support sigma in {6,7,8} and 12:
for each support R and each flat Z (closed set of restricted Hadamard columns containing the ones column) whose
orthogonal complement contains a full-support vector, with t = 12 - |Z| <= 5: record (R, Z, dim of complement,
a basis of the complement).  Output: typeflats_' + sys.argv[1].replace(',', '_') + '.pkl and a summary."""
import itertools, pickle
import numpy as np
from flint import fmpz_mat
from fplib import pal_patterns
M_, pats = pal_patterns(11)
core = [P for P in pats if len(P) == 6]
H = np.array([[1]*12] + [[1 if x in P else -1 for x in range(12)] for P in core], dtype=np.int64)
def rk(rows):
    return fmpz_mat(len(rows), len(rows[0]), [x for r in rows for x in r]).rank() if rows else 0
out = []
from collections import Counter
summ = Counter()
import sys
for s in [int(a) for a in sys.argv[1].split(",")]:
    for R in itertools.combinations(range(12), s):
        V = [[int(x) for x in H[c, list(R)]] for c in range(12)]
        seen = set(); stack = [[0]]
        while stack:
            Zl = stack.pop()
            r = rk([V[c] for c in Zl])
            if r >= s: continue
            Z = frozenset(c for c in range(12) if rk([V[x] for x in Zl] + [V[c]]) == r)
            if Z in seen: continue
            seen.add(Z)
            t = 12 - len(Z)
            if t <= int(sys.argv[2]):
                Mz = fmpz_mat(len(Z), s, [x for c in sorted(Z) for x in V[c]])
                X, nul = Mz.nullspace()
                if nul > 0 and all(any(X[i, j] != 0 for j in range(nul)) for i in range(s)):
                    basis = [[int(X[i, j]) for i in range(s)] for j in range(nul)]
                    out.append((s, t, R, tuple(sorted(Z)), nul, basis))
                    summ[(s, t, nul)] += 1
            if t > 5 - 0:  # deeper flats have smaller t: continue
                pass
            for c in range(12):
                if c not in Z: stack.append(sorted(Z) + [c])
    print('done support', s, flush=True)
for k, v in sorted(summ.items()):
    print('sigma %d t %d complement-dim %d : %d flats' % (k[0], k[1], k[2], v))
pickle.dump(out, open('typeflats_' + sys.argv[1].replace(',', '_') + '.pkl', 'wb'))
