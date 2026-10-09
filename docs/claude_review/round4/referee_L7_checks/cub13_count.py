"""NUMERICAL (referee): number of nondegenerate cub13 curves with x3 FIXED (z34, z35, p3 given, p36 from x3(1)=1).
Square system in 11 unknowns (z4 p45 p46 z5 p56 k5 c q345 q346 q356 q456) + barrier.  If the generic count
is 1, the family is rational over the x3-parameters (and each member over Q when they are rational)."""
import numpy as np, itertools, sys
from cub13lib import *
rng2 = np.random.default_rng(int(sys.argv[1]))
mode = sys.argv[2] if len(sys.argv) > 2 else 'rat'
if mode == 'rat':
    z34, z35, p3 = 2.0, -3.0, 5.0
else:
    z34, z35, p3 = rng2.normal(size=3) + 1j*rng2.normal(size=3)
p36 = 1 - (1-z34)*(1-z35)/(1-p3)
def full(y):
    return np.concatenate([[z34, z35, p3, p36], y])
def G(v):
    y = v[:11]; w = v[11:]; u = full(y)
    return np.concatenate([F(u)[1:], w*diffs(u) - 1])
found = []
for trial in range(int(sys.argv[3]) if len(sys.argv) > 3 else 300):
    y = rng2.normal(size=11)*2 + 1j*rng2.normal(size=11)*2
    v = np.concatenate([y, 1/diffs(full(y))])
    ok = False
    for it in range(80):
        g = G(v); ng = np.linalg.norm(g)
        if ng < 1e-12: ok = True; break
        dv = np.linalg.lstsq(jac(G, v), -g, rcond=None)[0]
        lam = 1.0
        while lam > 1e-4 and np.linalg.norm(G(v + lam*dv)) > ng: lam /= 2
        v = v + lam*dv
        if not np.all(np.isfinite(v)) or np.linalg.norm(v[:11]) > 1e5: break
    if not ok: continue
    u = full(v[:11])
    ev = collisions(u); cnt = {}; bad = 0
    for r, big in ev:
        if len(big) != 1: bad += 1; continue
        k = canon(big[0]); cnt[k] = cnt.get(k, 0) + 1
    good = (bad == 0 and set(cnt) == expected | {frozenset((1, 6))} and cnt.get(frozenset((1, 6)), 0) == 3
            and all(val == 1 for k, val in cnt.items() if k != frozenset((1, 6))))
    if not good: continue
    # identify up to the S_? symmetry: compare (z4,p45,p46,z5,p56,k5,c) as a set-invariant key
    key = np.round(np.array([u[4], u[9], u[10]]), 6)
    if not any(np.allclose(key, f[0], atol=1e-5) for f in found):
        found.append((key, u))
        print('new solution', len(found), 'at trial', trial, ':', np.round(u[4:11], 6), flush=True)
print('distinct nondegenerate solutions:', len(found))
np.save('cub13_count_%s_%s.npy' % (sys.argv[1], mode), np.array([f[1] for f in found]))
