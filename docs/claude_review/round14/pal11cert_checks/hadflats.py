"""Exact: for every support R (subset of the 12 rows, |R| >= 2) and every nonzero real z with supp z = R and
sum z = 0, the number t(z) of core columns c (11 Paley core columns) with <z, beta_c> != 0.
Computes f(R) = min t(z) over such z, by enumerating the subsets Z of the 12 Hadamard columns (all-ones included,
forced) that are 'closed' flats of the restricted vectors h_c|_R, and checking whether the orthogonal complement of
span{h_c|_R : c in Z} contains a vector of full support R (generic element of the complement has full support iff
the complement lies in no coordinate hyperplane).  Output: min t as a function of |R|, and the extremal supports."""
import itertools
import numpy as np
from sympy import Matrix
from fplib import pal_patterns
M_, pats = pal_patterns(11)
core = [P for P in pats if len(P) == 6]
# Hadamard columns as +-1 vectors on rows 0..11: beta_c = 2*1_P - 1 (P relative pattern, base row 0 not in P)
H = [[1]*12] + [[1 if x in P else -1 for x in range(12)] for P in core]
H = np.array(H, dtype=np.int64)   # 12 x 12 (rows = columns of the Hadamard matrix)
assert (H @ H.T == 12 * np.eye(12, dtype=np.int64)).all()
from fractions import Fraction
def rank(vs):
    return Matrix(vs).rank() if vs else 0
best = {}
ext = {}
for s in range(2, 13):
    mn = 99
    for R in itertools.combinations(range(12), s):
        V = H[:, list(R)]   # 12 x s: restricted vectors
        # enumerate flats greedily: a set Z containing column 0 (ones); its closure = all c with V[c] in span(V[Z]).
        # Max |Z| over flats with complement having a full-support vector.  Enumerate by DFS over closures.
        seen = set()
        def closure(Zl):
            Mz = Matrix(V[list(Zl)].tolist())
            r = Mz.rank()
            return frozenset(c for c in range(12) if Matrix(V[list(Zl) + [c]].tolist()).rank() == r), r
        stack = [closure([0])]
        localmax = -1
        while stack:
            Z, r = stack.pop()
            if Z in seen or r >= s: continue
            seen.add(Z)
            # complement: nullspace of V[Z] (as rows) in R^s
            ns = Matrix(V[sorted(Z)].tolist()).nullspace()
            if ns:
                full = all(any(v[i] != 0 for v in ns) for i in range(s))
                if full:
                    if len(Z) > localmax: localmax = len(Z)
            for c in range(12):
                if c not in Z:
                    stack.append(closure(sorted(Z) + [c]))
        if localmax >= 0:
            t = 12 - localmax   # core columns not in Z (Z contains the ones column)
            if t < mn: mn = t; ext[s] = [R]
            elif t == mn: ext[s].append(R)
    best[s] = mn
    print('support', s, ': min t_core =', mn, ' #extremal supports', len(ext.get(s, [])), flush=True)
