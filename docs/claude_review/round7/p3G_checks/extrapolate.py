"""Evidence/extrapolation: divisor density of L_X by a log-grid DP, calibrated against the measured
tight pair degree (threshold_out.jsonl), and the predicted closure threshold K*(M) for larger M.
Not a proof; the DP rounds logs to a grid of step h (counts are approximate)."""
import json, math, numpy as np
from gtight import primes_upto
h = 0.01
def log_divisor_counts(X, tmax):
    nb = int(tmax / h) + 2
    c = np.zeros(nb); c[0] = 1.0
    for q in primes_upto(X):
        if q == 2: continue
        old = c.copy()
        Q = q
        while Q <= X:
            sh = int(round(math.log(Q) / h))
            if sh < nb:
                c[sh:] += old[:nb - sh]
            Q *= q
    return c
def density(c, logY):
    # divisors in [Y/2, 2Y] per odd integer in that window
    lo, hi = int((logY - math.log(2)) / h), int((logY + math.log(2)) / h)
    return c[lo:hi + 1].sum() / (0.75 * math.exp(logY))
rows = [json.loads(l) for l in open('threshold_out.jsonl')]
# calibration: measured pair degree vs 8(M-1)*(6/8)*density(log(c P))
cache = {}
def dens(M, logY):
    X = 2 * M
    if X not in cache:
        cache[X] = log_divisor_counts(X, 8 * math.log(M) + 5)
    return density(cache[X], logY)
best = None
used = [d for d in rows if d['deg'] >= 0.5]
for lc in np.arange(-1.5, 1.01, 0.05):
    err = 0
    for d in used:
        M, K = d['M'], d['K']
        pred = 8 * (M - 1) * 0.75 * dens(M, K * math.log(M) + lc)
        err += (math.log(pred) - math.log(d['deg'])) ** 2
    if best is None or err < best[0]: best = (err, lc)
print('calibration: log c = %.2f, rms log error = %.3f over %d runs (pair degree >= 0.5)' % (best[1], math.sqrt(best[0] / len(used)), len(used)))
lc = best[1]
print('| M | predicted K* (pair degree 7.5) | lambda* = K* log M | lambda*/H | K*/(log M/loglog M) |')
print('|---|---|---|---|---|')
for M in [84, 164, 324, 644, 1284, 10**4, 10**5, 10**6]:
    if 2 * M not in cache:
        cache[2 * M] = log_divisor_counts(2 * M, 12 * math.log(M) + 5)
    lo, hi = 2.0, 12.0
    for _ in range(40):
        K = (lo + hi) / 2
        deg = 8 * (M - 1) * 0.75 * dens(M, K * math.log(M) + lc)
        if deg > 7.5: lo = K
        else: hi = K
    L = math.log(M); H = L * L / math.log(L)
    print(f'| {M} | {K:.2f} | {K*L:.1f} | {K*L/H:.2f} | {K/(L/math.log(L)):.2f} |')
