"""Equal-weight SA on genuine targets from a structured start: k stacked Paley copies plus random columns."""
import sys, math, time, numpy as np
sys.path.insert(0, '..')
from eqw import *
from gauss7 import paley_normalised
from qcore import primes_1mod4_above
M, t, Xf = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
D = np.load('../data/D_construct_M%d_t%d_u1_X%d.npy' % (M, t, Xf))
X = M * Xf; w0 = math.log(primes_1mod4_above(X, 1)[0]); lam2 = 2 * np.linalg.eigvalsh(D)[-1]
P = paley_normalised(M - 1)[:, 1:].astype(np.int8)
rng = np.random.default_rng(7)
for f in (1.0, 1.5, 2.0, 3.0):
    tau = f * w0; rP = plotkin_r(D, tau); found = None
    for r in range(rP, rP + 60, 2):
        T = thresholds(D, tau, r); k = r // (M - 1)
        best = None
        for j in range(2):
            perm = rng.permutation(M)
            S0 = np.hstack([P[perm]] * k + [random_balanced(M, r - k * (M - 1), rng)]) if k else random_balanced(M, r, rng)
            v, S = sa_run(T, S0, 10_000_000, T0=0.5, seed=j + 11)
            best = v if best is None else min(best, v)
            if v == 0 and exact_ok(S, T):
                found = (r, S); break
        if found: break
    if found:
        pol = polish(D, found[1], w0)[0]
        print('M=%d t=%d X=%d f=%.1f rP=%d r=%d r*tau=%.1f polish=%.1f ratio %.3f' % (M, t, X, f, rP, found[0], found[0]*tau, pol, pol/lam2), flush=True)
    else:
        print('M=%d t=%d X=%d f=%.1f rP=%d: none up to rP+58 (last V %s)' % (M, t, X, f, rP, best), flush=True)
