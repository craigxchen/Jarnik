"""NUMERICAL (referee): as cub13_frame.py, with barrier unknowns w_k and equations w_k * d_k(u) = 1 for the
differences d_k that collapse in the degenerate (cancellation) solutions: zero-pole pairs of each mover, c,
k5-1, and q's away from each other and from the special points."""
import numpy as np, itertools, sys
exec(open('cub13_frame.py').read().split('sols = []')[0])
def diffs(u):
    z34, z35, p3, p36, z4, p45, p46, z5, p56, k5, c, q345, q346, q356, q456 = u
    special = [z34, z35, p3, p36, z4, p45, p46, z5, p56, 0, 1]
    d = [z34-p3, z34-p36, z35-p3, z35-p36, z4-p45, z4-p46, z34-p45, z34-p46, z5-p45, z5-p56, z35-p45, z35-p56,
         c, k5-1, z34-z35, p3-p36, p45-p46, p45-p56, p46-p56, p36-p46, p36-p56, z4-z34, z5-z35, z4-z5]
    qs = [q345, q346, q356, q456]
    d += [a-b for a, b in itertools.combinations(qs, 2)]
    d += [qq - sp for qq in qs for sp in special]
    return np.array(d)
ND = len(diffs(np.ones(15) + np.arange(15)*0.1j))
def G(v):
    u = v[:15]; w = v[15:]
    return np.concatenate([F(u), w*diffs(u) - 1])
sols = []
for trial in range(400):
    u = rng.normal(size=15) + 1j*rng.normal(size=15)
    v = np.concatenate([u, 1/diffs(u)])
    ok = False
    for it in range(80):
        g = G(v); ng = np.linalg.norm(g)
        if ng < 1e-12: ok = True; break
        dv = np.linalg.lstsq(jac(G, v), -g, rcond=None)[0]
        lam = 1.0
        while lam > 1e-4 and np.linalg.norm(G(v + lam*dv)) > ng: lam /= 2
        v = v + lam*dv
        if not np.all(np.isfinite(v)) or np.linalg.norm(v[:15]) > 1e5: break
    if not ok: continue
    u = v[:15]
    ev = collisions(u)
    cnt = {}; bad = 0
    for r, big in ev:
        if len(big) != 1: bad += 1; continue
        k = canon(big[0]); cnt[k] = cnt.get(k, 0) + 1
    good = (bad == 0 and set(cnt) == expected | {frozenset((1, 6))} and cnt.get(frozenset((1, 6)), 0) == 3
            and all(val == 1 for k, val in cnt.items() if k != frozenset((1, 6))))
    Jm = jac(F, u); sv = np.linalg.svd(Jm, compute_uv=False)
    _, S_, Vh = np.linalg.svd(Jm); N = Vh[12:].conj().T
    s0 = 0.37 + 0.21j
    def ev_map(x):
        uu = u + N @ x[:3]; return np.array(X(uu, s0 + x[3]))
    se = np.linalg.svd(jac(ev_map, np.zeros(4, complex)), compute_uv=False)
    sols.append(u)
    print('trial %3d nondeg=%s events %d bad %d  rank-margin svJ %.2e  eval-map sv %.2e  events:%s'
          % (trial, good, len(ev), bad, sv[-1]/sv[0], se[-1]/se[0], sorted((tuple(sorted(k)), n) for k, n in cnt.items()) if not good else ''), flush=True)
    if len(sols) >= 6: break
print('converged', len(sols), 'of', trial+1)
np.save('cub13_frame_b_sols.npy', np.array(sols))
