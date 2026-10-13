import numpy as np, sys
from q4 import *
rng = np.random.default_rng(3)
conv = 0
for trial in range(60):
    u0 = rng.normal(size=18) + 1j * rng.normal(size=18)
    u = newton(u0, [0, 5, 8, 12])
    if u is None: continue
    conv += 1
    ok, cnt, bad, md = genuine(u)
    diff = {tuple(sorted(k)): v for k, v in cnt.items() if EXPECT.get(k, 0) != v}
    miss = {tuple(sorted(k)): v for k, v in EXPECT.items() if cnt.get(k, 0) != v}
    print('conv', ok, 'md %.2e' % md, 'nbad', len(bad), 'extra/diff', diff, 'expected-mismatch', miss)
    if conv > 6: break
print('converged', conv)
