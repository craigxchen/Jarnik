"""Sign structure of the 5 sdc11 tangent octics (evidence)."""
import json, numpy as np
D = json.load(open('sdc11_m8_basis.json'))
A = np.array(D['A'], float); C = np.array(D['basis'], float); C /= np.abs(C).sum(1, keepdims=True)
rng = np.random.default_rng(0)
U = rng.uniform(-1, 1, (4000, 11))
Y = (U @ A.T) ** 8
P = Y @ C.T
print('random u in [-1,1]^11: fraction positive per form:', np.round((P > 0).mean(0), 3))
# well separated: random permutations of equally spaced points
for s_frac in [1.0, 0.7, 0.4]:
    pts = []
    for _ in range(3000):
        base = np.linspace(-1, 1, 12)
        jitter = rng.uniform(-1, 1, 12) * (1 - s_frac) * (2 / 11) / 2
        z = base + jitter
        rng.shuffle(z)
        u = z[1:] - z[0]
        u /= np.abs(u).max()
        pts.append(u)
    U2 = np.array(pts)
    P2 = ((U2 @ A.T) ** 8) @ C.T
    print('near-equally spaced (jitter %.1f): fraction positive per form:' % (1 - s_frac), np.round((P2 > 0).mean(0), 3),
          ' all-same-sign rows: %.3f' % (np.all(P2 > 0, 1) | np.all(P2 < 0, 1)).mean())
