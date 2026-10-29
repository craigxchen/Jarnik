"""J9a' (rsyl:4:1,...,8,0:3): search ternary forms phi with |phi(A)| = R and |phi(A_l)| <= R-1 for all l in L_pair
(Thm 2.3 range condition, both signs of each shift), for R-1 = 9 (the failing degree). Also report sigma's range/rigidity."""
import sys, itertools
import numpy as np
sys.path.insert(0, '/tmp/claude-0/-home-user-Jarnik/ba4929a4-8665-5790-985a-3f09130a14b9/scratchpad/round28/verify_checks')
from vlib import *
spec = sys.argv[1] if len(sys.argv) > 1 else 'rsyl:4:1,2,3,4,5,6,7,8,0:3'
target = int(sys.argv[2]) if len(sys.argv) > 2 else 9
n, pats, L, A = spec_box(spec); d = n - 1
idx = {tuple(k.tolist()): i for i, k in enumerate(A)}
masks = []
for l in lpair_classes(d):
    for s in (1, -1):
        Al = shifted(A, np.array(l) * s, idx)
        m = np.zeros(len(A), bool); m[[idx[tuple(k.tolist())] for k in Al]] = True
        masks.append(m)
Ai = A.astype(np.int64)
hits = []
for phi in itertools.product([-1, 0, 1], repeat=d):
    if not any(phi): continue
    if next(x for x in phi if x) < 0: continue
    v = Ai @ np.array(phi)
    R = len(np.unique(v))
    if R - 1 != target: continue
    if all(len(np.unique(v[m])) <= R - 1 for m in masks):
        hits.append(phi)
sig = Ai.sum(axis=1)
print(spec, 'sigma range', sig.min(), sig.max(), '; forms phi with range', target + 1, 'on A and <=', target, 'on every A_l:', len(hits))
for h in hits[:20]: print('  ', h)
