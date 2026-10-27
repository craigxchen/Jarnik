"""H_light check for sdc_q: exact minimum of sigma_lambda / r over all integral sum-zero lambda with ||lambda||_1 = 4
that are NOT supported on two points (types e_x+e_y-e_z-e_w and 2e_x-e_y-e_z), sdc_q, mu = 1.  Exact (Fractions)."""
import itertools, time
from fractions import Fraction
from flib import sdc_patterns
import numpy as np
for q in (7, 11, 19, 23, 31):
    t0 = time.time()
    M, pats = sdc_patterns(q)
    r = len(pats)
    I = np.zeros((r, M), dtype=np.int64)
    for i, P in enumerate(pats):
        I[i, list(P)] = 1
    best4 = None; best2 = None
    # type e_x+e_y-e_z-e_w : vectorised over w
    for x, y in itertools.combinations(range(M), 2):
        for z in range(M):
            if z in (x, y): continue
            for w in range(z + 1, M):
                if w in (x, y): continue
                v = I[:, x] + I[:, y] - I[:, z] - I[:, w]
                s = Fraction(int(np.abs(v).sum() * 2 - r), 2)
                if best4 is None or s < best4[0]: best4 = (s, (x, y, z, w))
    for x in range(M):
        for y, z in itertools.combinations([t for t in range(M) if t != x], 2):
            v = 2 * I[:, x] - I[:, y] - I[:, z]
            s = Fraction(int(np.abs(v).sum() * 2 - r), 2)
            if best2 is None or s < best2[0]: best2 = (s, (x, y, z))
    print(f"sdc{q}: r={r}  min sigma, type ++--: {best4[0]} (= {float(best4[0])/r:.3f} r)  "
          f"type 2--: {best2[0]} (= {float(best2[0])/r:.3f} r)   {time.time()-t0:.0f}s", flush=True)
