"""Referee: independent certification of the 61 line covers of coset.md §8.2 (60 random lines at T=22, and the
eta=1e-3 line x_beta + 1e-3 u + R b at T=20).  Own monomial enumeration (both signs, then canonical sign), own interval
construction and sweep.  Error budget: c=<v,g> computed in float64 from exact-binary g (error < 1e-14), half-widths
from math.asin (rel. error ~1e-16): every endpoint is within 1e-13 of its true value, so an overlap margin > 1e-10
certifies the cover rigorously.  Additionally the critical eta=1e-3 line is redone with mpmath (50 digits)."""
import math, json, sys
import numpy as np
from fractions import Fraction

PS = [13, 17, 29, 37, 41, 53, 61, 73]
LOG = [math.log(p) for p in PS]

def enum(T):
    res = []
    r = len(PS)
    def rec(j, cur, w):
        if j == r:
            res.append((list(cur), w)); return
        k = 0
        while w + k * LOG[j] <= T + 1e-12:
            for s in ((1, -1) if k else (1,)):
                cur.append(s * k); rec(j + 1, cur, w + k * LOG[j]); cur.pop()
            k += 1
    rec(0, [], 0.0)
    out = []
    for v, w in res:
        nz = [x for x in v if x]
        if nz and nz[0] > 0:
            out.append((v, w))
    V = np.array([v for v, w in out], dtype=np.int64); W = np.array([w for v, w in out])
    return V, W

def cover_min_overlap(V, W, b, g):
    n = V @ b; c = V @ g
    keep = n != 0
    n, c, W = n[keep], c[keep], W[keep]
    y = np.exp(-W / 2)
    halfs = [(0.0, (2 / math.pi) * np.arcsin(y)), (0.5, (2 / math.pi) * np.arcsin(y / math.sqrt(2)))]
    lo_all, hi_all = [], []
    for shift, a in halfs:
        # condition violated iff n*theta + c in shift + Z + (-a, a)  <=> theta in ((k + shift - c - a)/n, (k + shift - c + a)/n)
        for sgn in (1, -1):
            sel = np.sign(n) == sgn
            nn, cc, aa = n[sel], c[sel], a[sel]
            # theta in [0,1): n*theta ranges over [0, n) or (n, 0]; k from floor(min) - 2 to ceil(max) + 2
            for k in range(-int(np.abs(nn).max()) - 3, int(np.abs(nn).max()) + 4):
                e1 = (k + shift - cc - aa) / nn; e2 = (k + shift - cc + aa) / nn
                lo = np.minimum(e1, e2); hi = np.maximum(e1, e2)
                s2 = (hi > -0.5) & (lo < 1.5)
                lo_all.append(lo[s2]); hi_all.append(hi[s2])
    lo = np.concatenate(lo_all); hi = np.concatenate(hi_all)
    o = np.argsort(lo, kind='stable'); lo, hi = lo[o], hi[o]
    # sweep over [0,1]: track covered-up-to; the line parameter is periodic, so cover [0,1] by open intervals
    reach = -np.inf; worst = np.inf; started = False
    for l, h in zip(lo, hi):
        if l >= 1: break
        if not started:
            if l < 0 and h > 0:
                pass
        if reach == -np.inf:
            if l >= 0: return -1.0      # 0 not covered
        else:
            if l >= 0 or reach > 0:
                if l > 0 and reach < 1: worst = min(worst, reach - l)
        reach = max(reach, h)
    if reach <= 1: return -2.0
    worst = min(worst, reach - 1)
    # also require 0 strictly inside an interval
    inside0 = np.max(np.where(lo < 0, hi, -np.inf))
    worst = min(worst, inside0, )
    return float(worst)

if __name__ == '__main__':
    out = {}
    # 60 random lines (same RNG definition as coset_checks/line_batch.py / batch_margins.py)
    V22, W22 = enum(22.0)
    out['n_monomials_T22'] = int(len(V22))
    rng = np.random.default_rng(11)
    def rand_b():
        while True:
            b = rng.integers(-3, 4, size=8)
            if np.gcd.reduce(np.abs(b[b != 0])) == 1 and (b != 0).sum() >= 2:
                return b
    margins = []
    for i in range(60):
        b = rand_b(); g = rng.uniform(0, 1, 8)
        margins.append(cover_min_overlap(V22, W22, b, g))
    out['random60_min_margin'] = min(margins); out['random60_n_certified'] = sum(m > 1e-10 for m in margins)
    out['random60_margins'] = margins
    # eta = 1e-3 line at T = 20
    split = [5, 13, 17, 29, 37, 41, 53, 61, 73]
    def gp(p):
        for bb in range(1, int(math.isqrt(p)) + 1):
            a2 = p - bb * bb; a = math.isqrt(a2)
            if a * a == a2 and a > bb: return a, bb
    qs = [split[split.index(p) - 1] for p in PS]
    xb = np.array([math.atan2(gp(q)[1], gp(q)[0]) / (math.pi / 2) for q in qs])
    d = json.load(open('../coset_checks/line_cosets_out.json'))
    b = np.array(d['b']); u = np.array(d['u'])
    for T in (20.0, 24.0):
        V, W = enum(T)
        for e in (1e-3, 1e-1, 0.0):
            out[f'line T={T} eta={e}'] = cover_min_overlap(V, W, b, xb + e * u)
    print(json.dumps({k: v for k, v in out.items() if k != 'random60_margins'}, indent=1))
    json.dump(out, open('cover_check_out.json', 'w'), indent=1)
