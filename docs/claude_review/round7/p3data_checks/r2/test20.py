import sys, math, time, numpy as np
sys.path.insert(0, '..')
from eqw import *
from qcore import primes_1mod4_above
rng = np.random.default_rng(2)
M = int(sys.argv[1]); ts = [int(v) for v in sys.argv[2].split(',')]; Xf = int(sys.argv[3]); iters = int(float(sys.argv[4]))
for t in ts:
    D = np.load('../data/D_construct_M%d_t%d_u1_X%d.npy' % (M, t, Xf))
    X = M * Xf; w0 = math.log(primes_1mod4_above(X, 1)[0])
    lam = np.linalg.eigvalsh(D)[-1]
    t0 = time.time(); best = None
    for f in (1.0, 1.25, 1.5, 1.75, 2.0, 2.5, 3.0):
        tau = f * w0
        rP, r, S = search_tau(D, tau, rng, iters=iters, restarts=2, extra=30)
        pol = polish(D, S, w0)[0] if S is not None else None
        if pol is not None and (best is None or pol < best[0]): best = (round(pol, 2), f, r, rP)
        print('  t=%d f=%.2f rP=%d (rP tau %.1f) r=%s rtau=%s polish=%s  %.0fs' % (t, f, rP, rP*tau, r, None if r is None else round(r*tau,2), None if pol is None else round(pol,2), time.time()-t0), flush=True)
    print('M=%d X=%d t=%d 2lam=%.2f best %s ratio %.3f (%.0fs)' % (M, X, t, 2*lam, best, best[0]/(2*lam), time.time()-t0), flush=True)
