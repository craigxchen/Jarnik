"""NUMERICAL (referee): collect the distinct NONDEGENERATE cub13 curves through the fixed point X0 found by
monodromy (both runs), re-verify each: residual of the 16x16 system, smallest singular value of its
Jacobian, 19 distinct contact points each lying on exactly one predicted divisor (D_16 three times)."""
import numpy as np, itertools
from cub13lib import *
u0 = np.load('cub13_frame_b_sols.npy')[0]; s0 = 0.31 + 0.47j; X0 = np.array(X(u0, s0))
def H(v):
    return np.concatenate([F(v[:15]), np.array(X(v[:15], v[15])) - X0])
cands = list(np.load('cub13_monodromy_sols_seed1.npy')) + list(np.load('cub13_monodromy_sols_s7.npy'))
good = []
for v in cands:
    # Newton polish
    for _ in range(5):
        v = v - np.linalg.lstsq(jac(H, v), H(v), rcond=None)[0]
    u = v[:15]
    z34, z35, p3, p36, z4, p45, p46, z5, p56, k5, c, q345, q346, q356, q456 = u
    num, den = nd(u)
    z6 = np.roots(np.trim_zeros(num[6], 'f'))
    special = [z34, z35, p3, p36, z4, p45, p46, z5, p56, 0, 1, q345, q346, q356, q456] + list(z6)
    dmin = min(abs(a-b) for a, b in itertools.combinations(special, 2))
    ev = collisions(u, tol=1e-8); cnt = {}; bad = 0
    for r, big in ev:
        if len(big) != 1: bad += 1; continue
        k = canon(big[0]); cnt[k] = cnt.get(k, 0) + 1
    ok = (bad == 0 and len(ev) == 19 and set(cnt) == expected | {frozenset((1, 6))}
          and cnt[frozenset((1, 6))] == 3 and dmin > 1e-4)
    sv = np.linalg.svd(jac(H, v), compute_uv=False)
    if ok and not any(np.linalg.norm(v - w) < 1e-6 for w in good):
        good.append(v)
        print('curve %d: |H| %.1e  sv_min %.2e  sv_max %.2e  min special-point distance %.2e  s=%s'
              % (len(good), np.linalg.norm(H(v)), sv[-1], sv[0], dmin, np.round(v[15], 4)))
print('distinct nondegenerate cub13 curves through X0:', len(good))
np.save('cub13_through_X0.npy', np.array(good))
