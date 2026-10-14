"""Weight vectors (nonnegative integers on the tau<8 test directions) for the knapsack bounds of csearch3."""
import numpy as np
from fplib import pal_patterns
dat = np.load('dirs_pal11.npz'); D = dat['D']; tau = dat['tau']
M_, pats = pal_patterns(11)
dirindex = {tuple(v): j for j, v in enumerate(D.tolist())}
def ind(P):
    v = [0]*11
    for x in P: v[x-1] = 1
    return dirindex[tuple(v)]
core = [ind(P) for P in pats if len(P) == 6]
twin = [ind(P) for P in pats if len(P) == 5]
coords = [j for j in range(len(D)) if tau[j] == 2]
ones = dirindex[tuple([1]*11)]
assert len(core) == 11 and len(twin) == 12 and len(coords) == 11
vecs = []
def vec(js, w=1):
    y = [0]*len(D)
    for j in js: y[j] = w
    return y
vecs.append(vec(core)); vecs.append(vec(twin)); vecs.append(vec(core + twin))
vecs.append(vec(coords + [ones]))
vecs.append(vec([j for j in range(len(D)) if tau[j] == 4]))
vecs.append(vec([j for j in range(len(D)) if tau[j] == 6]))
vecs.append([1]*len(D))
# mixtures: core*k + coords
for k in (1, 2, 3):
    y = vec(core, k)
    for j in coords + [ones]: y[j] += 1
    vecs.append(y)
with open('weights.txt', 'w') as f:
    f.write(f'{len(vecs)}\n')
    for y in vecs: f.write(' '.join(map(str, y)) + '\n')
print('core', core, 'twin', twin, 'ones', ones, 'nvecs', len(vecs))
