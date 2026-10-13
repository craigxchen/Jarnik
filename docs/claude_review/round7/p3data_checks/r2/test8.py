import sys, math, time, numpy as np
sys.path.insert(0, '..')
from eqw import *
from qcore import primes_1mod4_above
rng = np.random.default_rng(1)
for Xf in (2, 4):
    for t in range(6):
        D = np.load('../data/D_construct_M8_t%d_u1_X%d.npy' % (t, Xf))
        X = 8 * Xf; w0 = math.log(primes_1mod4_above(X, 1)[0])
        t0 = time.time(); best = None
        for f in (1.0, 1.25, 1.5, 2.0, 2.5, 3.0):
            tau = f * w0
            rP, r, S = search_tau(D, tau, rng, iters=1_000_000, restarts=2, extra=20)
            pol = polish(D, S, w0)[0] if S is not None else None
            if pol is not None and (best is None or pol < best[0]): best = (pol, f, r, rP)
            print('  X=%d t=%d f=%.2f rP=%d r=%s rtau=%s polish=%s' % (X, t, f, rP, r, None if r is None else round(r*tau,2), None if pol is None else round(pol,2)))
        print('X=%d t=%d best %s  (%.0fs)' % (X, t, best, time.time()-t0), flush=True)
