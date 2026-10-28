"""Robustness margins for the line covers of line_cosets.py, and location of survivors on the eta=0 line.
A cover is certified if consecutive sorted open intervals overlap by > 1e-9 everywhere on [0,1] (float error ~1e-13)."""
import sys, math, json
import numpy as np
sys.path.insert(0, '.')
from common import *
from line_cosets import *

def intervals(V, W, b, g):
    n = V @ b; c = V @ g; a0, a1 = halfwidths(W)
    nc = n != 0; n = n[nc]; c = c[nc]; a0 = a0[nc]; a1 = a1[nc]
    sgn = np.sign(n); n = np.abs(n); c = c * sgn
    los, his = [], []
    for shift, a in ((0.0, a0), (0.5, a1)):
        cc = (c - shift) % 1.0
        for k in range(-1, n.max() + 3):
            sel = (k <= n + 1)
            los.append((k - cc[sel] - a[sel]) / n[sel]); his.append((k - cc[sel] + a[sel]) / n[sel])
    lo = np.concatenate(los); hi = np.concatenate(his)
    keep = (hi > -0.01) & (lo < 1.01)
    lo = lo[keep]; hi = hi[keep]; o = np.argsort(lo)
    return lo[o], hi[o]

def cover_margin(lo, hi):
    # need: some interval contains 0 with margin, consecutive overlaps, and reach past 1
    hmax = np.maximum.accumulate(hi)
    first = lo[0]
    over = hmax[:-1] - lo[1:]           # >0 means overlap (only relevant where lo[1:] within [0,1])
    rel = (lo[1:] > 0) & (lo[1:] < 1)
    return float(min(-first, over[rel].min(), hmax[-1] - 1))

allsp = split_primes(200)
ps = [13, 17, 29, 37, 41, 53, 61, 73]
qs = [allsp[allsp.index(p) - 1] for p in ps]
xb = np.array([math.atan2(*gaussian_prime(q)[::-1]) / (math.pi / 2) for q in qs])
logs = [math.log(p) for p in ps]
d = json.load(open('line_cosets_out.json'))
b = np.array(d['b']); u = np.array(d['u'])
out = {}
for T in (20.0, 24.0):
    V, W = monomials(logs, T)
    for e in (1e-3, 1e-1):
        lo, hi = intervals(V, W, b, xb + e * u)
        out[f'T={T},eta={e}'] = cover_margin(lo, hi)
    lo, hi = intervals(V, W, b, xb)
    # survivors on eta=0 line: report gap locations theta (x = xb + theta b), theta in [0,1)
    hmax = np.maximum.accumulate(hi); prev = np.concatenate(([lo[0]], hmax[:-1]))
    gi = np.nonzero((lo > prev) & (lo > 0) & (prev < 1))[0]
    out[f'T={T},eta=0,gaps'] = [(float(prev[i]), float(lo[i])) for i in gi]
print(json.dumps(out, indent=1))
json.dump(out, open('line_verify_out.json', 'w'), indent=1)
