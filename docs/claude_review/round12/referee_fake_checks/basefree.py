import numpy as np
from rlib import *
for name, n, pats, L in [('sdc7', 8, sdc(7), 2), ('sdc11', 12, sdc(11), 2), ('pal7', 8, pal(7), 2), ('pal7', 8, pal(7), 3), ('pal11', 12, pal(11), 2)]:
    A, others = box(n, pats, L, 0)
    kb = -A.sum(1)
    sets = []
    for b in range(n):
        Ab, ob = box(n, pats, L, b)
        full = np.zeros((len(Ab), n), dtype=np.int64); full[:, ob] = Ab; full[:, b] = -Ab.sum(1)
        sets.append(set(map(tuple, full.tolist())))
    print(name, 'L', L, 'max |k_0| on A_0:', np.abs(kb).max(), ' all full-coordinate boxes equal:', all(s == sets[0] for s in sets), [len(s) for s in sets])
