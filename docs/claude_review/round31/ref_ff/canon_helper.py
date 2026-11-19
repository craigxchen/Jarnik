import itertools
_P=list(itertools.permutations(range(8)))
def canon_signs(mat):
    rows=[sum(1<<j for j in range(8) if r[j]==1) for r in mat]
    return canon_rows(rows)
def canon_rows(rows):
    import numpy as np
    best=None
    for fl in (0,255):
        R=[r^fl for r in rows]
        for p in _P:
            img=tuple(sorted(sum(((r>>j)&1)<<p[j] for j in range(8)) for r in R))
            if best is None or img<best: best=img
    return best
