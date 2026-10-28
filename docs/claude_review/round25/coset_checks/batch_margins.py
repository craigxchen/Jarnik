"""Cover margins (min overlap of consecutive open intervals) for the 60 random lines of line_batch.py (same RNG
sequence, T = 22). A positive margin >> 1e-12 certifies that the non-character windows cover the whole line."""
import sys, math, json
import numpy as np
sys.path.insert(0, '.')
from common import *
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
    keep = (hi > -0.01) & (lo < 1.01); lo = lo[keep]; hi = hi[keep]; o = np.argsort(lo)
    return lo[o], hi[o]
def cover_margin(lo, hi):
    hmax = np.maximum.accumulate(hi)
    over = hmax[:-1] - lo[1:]; rel = (lo[1:] > 0) & (lo[1:] < 1)
    return float(min(-lo[0], over[rel].min(), hmax[-1] - 1))
ps = [13, 17, 29, 37, 41, 53, 61, 73]
V, W = monomials([math.log(p) for p in ps], 22.0)
rng = np.random.default_rng(11)
def rand_b():
    while True:
        b = rng.integers(-3, 4, size=8)
        if np.gcd.reduce(np.abs(b[b != 0])) == 1 and (b != 0).sum() >= 2:
            return b
margins = []
for i in range(60):
    b = rand_b(); g = rng.uniform(0, 1, 8)
    lo, hi = intervals(V, W, b, g)
    margins.append(cover_margin(lo, hi))
print('min margin', min(margins), 'max', max(margins), 'n positive', sum(m > 1e-9 for m in margins))
json.dump(margins, open('batch_margins_out.json', 'w'))
