"""NUMERICAL (referee): how many cub13 curves pass through a general point of M_{0,7}?
Start: the nondegenerate solution u0 (cub13_frame_b_sols.npy) and a point s0 on it; target X0 = X(u0, s0).
System in (u, s): F(u) = 0 (12 eqs) and X(u, s) = Xtarget (4 eqs): 16 x 16.
Monodromy: move Xtarget around random triangle loops X0 -> X1 -> X2 -> X0 (complex), track the solution
by predictor-corrector, and record every distinct solution reached at X0.  If many loops never produce a
second solution, the count is (numerically) 1, i.e. the curve through a general point is unique, hence
defined over Q when the point is."""
import numpy as np, sys
from cub13lib import F, X, jac, collisions, canon, expected
rng = np.random.default_rng(int(sys.argv[1]) if len(sys.argv) > 1 else 0)
u0 = np.load('cub13_frame_b_sols.npy')[0]
s0 = 0.31 + 0.47j
X0 = np.array(X(u0, s0))

def H(v, Xt):
    u, s = v[:15], v[15]
    return np.concatenate([F(u), np.array(X(u, s)) - Xt])

def newton(v, Xt, it=8, tol=1e-12):
    for _ in range(it):
        h = H(v, Xt)
        if np.linalg.norm(h) < tol: return v, True
        v = v - np.linalg.lstsq(jac(lambda w: H(w, Xt), v), h, rcond=None)[0]
    return v, np.linalg.norm(H(v, Xt)) < 1e-9

def track(v, Xa, Xb, steps=400):
    dt = 1.0/steps; t = 0.0
    while t < 1 - 1e-12:
        h = min(dt, 1 - t)
        Xn = Xa + (t + h)*(Xb - Xa)
        # predictor: Euler step using dv/dt = -J^{-1} dH/dt, dH/dt = -(Xb - Xa) on last 4 eqs
        Jv = jac(lambda w: H(w, Xa + t*(Xb - Xa)), v)
        rhs = np.concatenate([np.zeros(12), (Xb - Xa)])
        vp = v + h*np.linalg.lstsq(Jv, rhs, rcond=None)[0]
        vn, ok = newton(vp, Xn)
        if ok and np.linalg.norm(vn - v) < 0.2*(1 + np.linalg.norm(v)):
            v = vn; t += h; dt = min(dt*1.5, 0.05)
        else:
            dt /= 2
            if dt < 1e-7: return None
    return v

sols = [np.concatenate([u0, [s0]])]
v_start = sols[0]
print('start residual', np.linalg.norm(H(v_start, X0)), 'sv', np.linalg.svd(jac(lambda w: H(w, X0), v_start), compute_uv=False)[[0,-1]])
nloops = int(sys.argv[2]) if len(sys.argv) > 2 else 20
fails = 0
for L in range(nloops):
    X1 = X0 + (rng.normal(size=4) + 1j*rng.normal(size=4))*0.8
    X2 = X0 + (rng.normal(size=4) + 1j*rng.normal(size=4))*0.8
    for base in list(sols):
        v = track(base, X0, X1)
        if v is None: fails += 1; continue
        v = track(v, X1, X2)
        if v is None: fails += 1; continue
        v = track(v, X2, X0)
        if v is None: fails += 1; continue
        if not any(np.linalg.norm(v - w) < 1e-6*(1 + np.linalg.norm(w)) for w in sols):
            sols.append(v)
            ev = collisions(v[:15]); cnt = {}; bad = 0
            for r, big in ev:
                if len(big) != 1: bad += 1; continue
                k = canon(big[0]); cnt[k] = cnt.get(k, 0) + 1
            good = (bad == 0 and set(cnt) == expected | {frozenset((1, 6))})
            print('loop', L, 'NEW solution #%d' % len(sols), 'nondegenerate:', good, flush=True)
    print('loop', L, 'solutions so far', len(sols), 'tracking failures', fails, flush=True)
    np.save('cub13_monodromy_sols_s%s.npy' % (sys.argv[1] if len(sys.argv) > 1 else '0'), np.array(sols))
print('final count of solutions through the point:', len(sols))
np.save('cub13_monodromy_sols.npy', np.array(sols))
