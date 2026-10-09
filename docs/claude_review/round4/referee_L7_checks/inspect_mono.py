import numpy as np, itertools, sys
from cub13lib import *
S = np.load(sys.argv[1])
for n, v in enumerate(S):
    u, s = v[:15], v[15]
    z34, z35, p3, p36, z4, p45, p46, z5, p56, k5, c, q345, q346, q356, q456 = u
    num, den = nd(u)
    z6 = np.roots(np.trim_zeros(num[6], 'f'))
    special = [z34, z35, p3, p36, z4, p45, p46, z5, p56, 0, 1, q345, q346, q356, q456] + list(z6)
    dmin = min(abs(a-b) for a, b in itertools.combinations(special, 2))
    ev = collisions(u, tol=1e-7)
    cnt = {}; bad = []
    for r, big in ev:
        if len(big) != 1: bad.append((r, big)); continue
        k = canon(big[0]); cnt[k] = cnt.get(k, 0) + 1
    print(n, 'max|special| %.2e  min dist %.2e  |c| %.2e  events %d  bad %s' % (max(abs(np.array(special))), dmin, abs(c), len(ev), [(np.round(r,3), b) for r, b in bad][:2]))
