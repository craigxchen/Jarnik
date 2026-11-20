#!/usr/bin/env python3
"""Random-angle model for F_delta(P).

Model: for each signed exponent vector e (b = prod pi^{e}), arg b mod pi/4 is uniform.
Expected count with Y0 <= h <= B:  sum_e  min(1, (8/pi) arcsin(min(1, h^{-delta/2}/2))).
DP histogram (bin width dt) of the number of signed vectors by log h, exact recurrence
y[t] = hist[t] + hist[t-L] + y[t-L] for the factor (1+x^L)/(1-x^L).

For large s we use the Rankin/saddle upper bound
log E <= log(4/pi) + min_{0<eta<sigma} [ -eta*T0 + sum_p log((1+p^{-(sigma-eta)})/(1-p^{-(sigma-eta)})) ].
"""
import math, sys
import numpy as np
from sympy import primerange

def hist_dp(ps, Tmax, dt):
    n = int(Tmax/dt) + 1
    h = np.zeros(n); h[0] = 1.0
    for p in ps:
        L = int(round(math.log(p)/dt))
        y = h.copy()
        # y[t] = h[t] + h[t-L] + y[t-L]
        y[L:] += h[:-L]
        for start in range(L, n, L):
            end = min(n, start+L)
            y[start:end] += y[start-L:end-L]
        h = y
    return h

def expected(ps, delta, gamma, logB, dt=0.01, logY0=None):
    if logY0 is None:
        logY0 = gamma*sum(math.log(p) for p in ps)
    h = hist_dp(ps, logB, dt)
    t = np.arange(len(h))*dt
    w = np.arcsin(np.minimum(1.0, 0.5*np.exp(-delta/2*t)))
    prob = np.minimum(1.0, w*8/math.pi)
    m = (t >= logY0) & (t <= logB)
    return float((h*prob)[m].sum())

def rankin(ps, delta, gamma, logY0=None):
    sig = delta/2
    lps = np.log(np.array(ps, dtype=float))
    T0 = gamma*lps.sum() if logY0 is None else logY0
    best = None
    for eta in np.linspace(1e-4, sig*(1-1e-4), 2000):
        x = np.exp(-(sig-eta)*lps)
        v = -eta*T0 + np.log((1+x)/(1-x)).sum()
        if best is None or v < best[0]: best = (v, eta)
    return best[0] + math.log(4/math.pi), best[1]

if __name__ == '__main__':
    mode = sys.argv[1]
    allp = [p for p in primerange(5, 3*10**7) if p % 4 == 1]
    if mode == 'small':
        for delta in (0.3, 0.4, 0.45, 0.49):
            for gamma in (0.1, 0.25, 0.4, 0.49):
                row = []
                for s in range(3, 13):
                    row.append('%.3g' % expected(allp[:s], delta, gamma, 2500.0, dt=0.005))
                print('delta=%.2f gamma=%.2f  E[F] s=3..12:' % (delta, gamma), ' '.join(row))
    elif mode == 'large':
        for delta in (0.3, 0.45, 0.49):
            for gamma in (0.1, 0.25, 0.49):
                out = []
                for s in (10, 30, 100, 300, 1000, 3000, 10**4, 3*10**4, 10**5, 3*10**5, 10**6):
                    v, eta = rankin(allp[:s], delta, gamma)
                    out.append('%d:%.3g' % (s, v/s))
                print('delta=%.2f gamma=%.2f  (log E)/s  ' % (delta, gamma), ' '.join(out))
    elif mode == 'sparse':
        # s primes = 1 mod 4 starting at Q
        for delta in (0.3, 0.49):
            for gamma in (0.1, 0.25):
                for Q in (10**3, 10**5, 10**7):
                    start = next(i for i, p in enumerate(allp) if p > Q)
                    out = []
                    for s in (3, 10, 30, 100, 300):
                        v, eta = rankin(allp[start:start+s], delta, gamma)
                        out.append('%d:%.3g' % (s, v))
                    print('delta=%.2f gamma=%.2f Q=%g  log E (Rankin bound):' % (delta, gamma, Q), ' '.join(out))
