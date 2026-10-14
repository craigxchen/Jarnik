"""Which collision hyperplanes divide the sdc11 octics / lie in TC?  (numerical test on random points)"""
import json, itertools, numpy as np
D = json.load(open('sdc11_m8_basis.json'))
A = np.array(D['A'], float); C = np.array(D['basis'], float); C /= np.abs(C).sum(1, keepdims=True)
rng = np.random.default_rng(2)
def P(U): return ((U @ A.T) ** 8) @ C.T
U = rng.uniform(-1, 1, (200, 11)); scale = np.abs(P(U)).mean(0)
print('typical |Phi_i|:', scale)
hyp = [(x, None) for x in range(11)] + list(itertools.combinations(range(11), 2))
for x, y in hyp:
    U = rng.uniform(-1, 1, (200, 11))
    if y is None: U[:, x] = 0
    else: U[:, y] = U[:, x]
    r = np.abs(P(U)).max(0) / scale
    tag = 'u_%d=0' % x if y is None else 'u_%d=u_%d' % (x, y)
    if (r < 1e-8).any():
        print(tag, 'forms vanishing identically:', np.where(r < 1e-8)[0].tolist(), ' max ratio others %.2e' % r[r >= 1e-8].min() if (r>=1e-8).any() else '')
# rank of the restricted system on each hyperplane
