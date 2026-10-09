"""Cross-check: sum_lam f^lam h_lam(d) (isotypic, from the sweep data) == direct Hilbert function."""
import json, sys
from combi import *
from direct import direct_profile
from modla import P1
data = {}
for M in (6, 7, 8):
    try:
        for line in open(f'data/shell_M{M}.jsonl'):
            r = json.loads(line)
            if r['prime'] == P1: data[(M, r['p'], r['n'])] = r
    except FileNotFoundError: pass
todo = [(6, t) for t in range(1, 9)] + [(7, t) for t in range(1, 8)] + [(8, t) for t in range(1, 7)]
for M, t in todo:
    for s in range(t % 2, t + 1, 2):
        p, n = (t + s) // 2, (t - s) // 2
        r = data.get((M, p, n))
        if r is None: continue
        dp = direct_profile(shell_points(M, p, n), P1)
        L = max(len(dp), max(len(c['profile']) for c in r['comps']))
        tot = [0] * L
        for c in r['comps']:
            f = hook_dim(tuple(c['lam']))
            pr = c['profile']
            for d in range(L): tot[d] += f * pr[min(d, len(pr) - 1)]
        dpp = [dp[min(d, len(dp) - 1)] for d in range(L)]
        print(M, s, t, '|S|', r['size'], 'direct reg', len(dp) - 1, 'iso reg', r['reg'], 'profiles equal', tot == dpp, flush=True)
