"""Compact Part-1 tables for p3data.md (from genuine_trials.jsonl and paley_selfcheck.json)."""
import json, math
import numpy as np
from collections import defaultdict
from demands7 import expected_mean_independent
recs = [json.loads(l) for l in open('genuine_trials.jsonl')]
g = defaultdict(list)
for r in recs:
    g[(r['M'], r['variant'], r['mode'])].append(r)
print('| M | X | E D (exact, indep.) | mean D | D_max median [5%,95%] | D_max/L | D_max/H | lambda/(M L) | max_x D_x/(M L) | (M-1)D_max/lambda |')
print('|---|---|---|---|---|---|---|---|---|---|')
for M in sorted({r['M'] for r in recs}):
    L = math.log(M); H = L * L / math.log(L)
    for Xf in (2, 4):
        st = [r['u1_X%d' % Xf] for r in g[(M, 'construct', 'window')] + g[(M, 'lower', 'window')] + g[(M, 'construct', 'iid')]]
        mx = np.array([s['max'] for s in st]); lam = np.array([s['lam'] for s in st]); mr = np.array([s['maxrow'] for s in st])
        print('| %d | %d | %.3f | %.3f | %.1f [%.1f, %.1f] | %.2f | %.2f | %.3f | %.3f | %.2f |' % (
            M, Xf * M, expected_mean_independent(M, Xf * M), np.mean([s['mean'] for s in st]), np.median(mx),
            np.percentile(mx, 5), np.percentile(mx, 95), np.median(mx) / L, np.median(mx) / H, lam.mean() / (M * L),
            mr.mean() / (M * L), np.median((M - 1) * mx / lam)))
print()
print('any-unit demands (construct window): M | X | E D_any (exact) | mean | D_max median | lambda/(M L)')
for M in sorted({r['M'] for r in recs}):
    L = math.log(M)
    for Xf in (2, 4):
        st = [r['any_X%d' % Xf] for r in g[(M, 'construct', 'window')]]
        print(M, Xf * M, '%.2f' % expected_mean_independent(M, Xf * M, anyunit=True), '%.2f' % np.mean([s['mean'] for s in st]),
              '%.1f' % np.median([s['max'] for s in st]), '%.2f' % np.mean([s['lam'] for s in st]) if False else '%.2f' % (np.mean([s['lam'] for s in st]) / (M * L)))
print()
print('number of trials per M: construct-window %d, lower-window %d, construct-iid %d' % (
    len(g[(12, 'construct', 'window')]), len(g[(12, 'lower', 'window')]), len(g[(12, 'construct', 'iid')])))
