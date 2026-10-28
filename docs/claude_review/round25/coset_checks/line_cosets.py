"""Section 8.2: one-dimensional cosets (lines) g + R b in T^r, r = 8 primes 13..73 (no 5).
For each line: exact union of all non-character windows with w(v) <= T (intervals in x in [0,1)),
character check (n = <b,v> = 0) at g, surviving length and number of surviving components.
g = transplant point xb + eta*u (u fixed random direction), and uniform random g."""
import sys, math, json
import numpy as np
sys.path.insert(0, '.')
from common import *

def line_survivors(V, W, b, g):
    n = V @ b
    c = V @ g
    a0, a1 = halfwidths(W)
    ch = (n == 0)
    char_margin = None
    if ch.any():
        m = np.minimum(dist_Z(c[ch]) - a0[ch], dist_Z(c[ch] - 0.5) - a1[ch])
        char_margin = float(m.min())
    nc = ~ch
    n = n[nc]; c = c[nc]; a0 = a0[nc]; a1 = a1[nc]
    sgn = np.sign(n); n = np.abs(n); c = c * sgn
    los, his = [], []
    for shift, a in ((0.0, a0), (0.5, a1)):
        cc = (c - shift) % 1.0            # condition: dist(n x + cc, Z) >= a
        kmax = n.max() + 2
        for k in range(-1, kmax + 1):
            sel = (k <= n + 1)
            lo = (k - cc[sel] - a[sel]) / n[sel]; hi = (k - cc[sel] + a[sel]) / n[sel]
            los.append(lo); his.append(hi)
    lo = np.concatenate(los); hi = np.concatenate(his)
    keep = (hi > 0) & (lo < 1)
    lo = np.clip(lo[keep], 0, 1); hi = np.clip(hi[keep], 0, 1)
    o = np.argsort(lo); lo = lo[o]; hi = hi[o]
    # sweep
    gaps = []
    cur = 0.0
    hmax = np.maximum.accumulate(hi)
    starts = lo
    prev = np.concatenate(([0.0], hmax[:-1]))
    gi = np.nonzero(starts > prev)[0]
    for i in gi:
        gaps.append((float(prev[i]), float(starts[i])))
    if hmax[-1] < 1:
        gaps.append((float(hmax[-1]), 1.0))
    # windows are open: a gap of positive length survives; length-0 touching ignored
    L = sum(h - l for l, h in gaps)
    return {'char_margin': char_margin, 'n_intervals': int(len(lo)), 'surv_len': L, 'n_gaps': len(gaps),
            'gaps_sample': gaps[:5], 'S_nc': float(window_measure(W[nc]).sum())}

if __name__ == '__main__':
    allsp = split_primes(200)
    ps = [13, 17, 29, 37, 41, 53, 61, 73]
    qs = [allsp[allsp.index(p) - 1] for p in ps]
    xb = np.array([math.atan2(*gaussian_prime(q)[::-1]) / (math.pi / 2) for q in qs])
    logs = [math.log(p) for p in ps]
    rng = np.random.default_rng(7)
    b = np.array([1, -1, 2, 0, 1, -2, 1, 1]); u = rng.normal(size=8); u /= np.linalg.norm(u)
    results = []
    for T in [16.0, 20.0, 24.0, 28.0]:
        V, W = monomials(logs, T)
        cases = [('eta=0', xb)] + [(f'eta={e:g}', xb + e * u) for e in (1e-9, 1e-7, 1e-5, 1e-3, 1e-1)]
        cases += [(f'uniform{i}', rng.uniform(0, 1, 8)) for i in range(4)]
        for name, g in cases:
            r = line_survivors(V, W, b, g)
            r.update({'T': T, 'case': name, 'n_mon': int(len(V))})
            results.append(r)
            print(json.dumps({k: r[k] for k in ('T', 'case', 'char_margin', 'surv_len', 'n_gaps', 'S_nc', 'n_intervals')}), flush=True)
    json.dump({'primes': ps, 'q': qs, 'b': b.tolist(), 'u': u.tolist(), 'results': results},
              open('line_cosets_out.json', 'w'), indent=1)
