"""Referee: exact re-check of M=8 actual-prime profiles (q_results_M8.jsonl): pairs, and for +sep the
Lemma 1.1(2) placement with rational positions; also W, prime properties; W_prime/W_sc ratios."""
import json, math
from fractions import Fraction
import numpy as np
from sympy import isprime
P8 = '/tmp/claude-0/-home-user-Jarnik/ba4929a4-8665-5790-985a-3f09130a14b9/scratchpad/round7/p3data_checks/q_results_M8.jsonl'
DATA = '/tmp/claude-0/-home-user-Jarnik/ba4929a4-8665-5790-985a-3f09130a14b9/scratchpad/round7/p3data_checks/data'
M = 8
def dec(mask): return [-1 if (mask >> x) & 1 else 1 for x in range(M)]
def check(cols, X, mod, t=None):
    S = [dec(m) for m, _ in cols]; P = [math.prod(pl) for _, pl in cols]
    allp = [p for _, pl in cols for p in pl]
    assert len(set(allp)) == len(allp) and all(isprime(p) and p % 4 == 1 and p > X for p in allp)
    assert all(abs(sum(s)) < M for s in S)
    rows = {tuple(s[x] for s in S) for x in range(M)}; assert len(rows) == M
    N = math.prod(P); ok = True; okp = True
    for x in range(M):
        for y in range(x+1, M):
            sep = math.prod(P[j] for j in range(len(P)) if S[j][x] != S[j][y]); same = N // sep
            m = int(mod[x, y])
            if not sep >= m**4*same: ok = False
            if t is not None:
                dt = abs(t[x] - t[y])
                if not dt**4*(sep - m*m) >= m**4*same: okp = False
    return ok, okp, sum(math.log(p) for p in allp), len(allp)
rat = []
for line in open(P8):
    R = json.loads(line); X, tt = R['X'], R['t']
    mod = np.load(f'{DATA}/mod_construct_M8_t{tt}_u1_X{X//8}.npy', allow_pickle=True)
    ok1, _, W1, r1 = check(R['prime_cols'], X, mod)
    t = [Fraction(a, b) for a, b in R['positions']]
    assert all(0 <= v <= 1 for v in t)
    ok2, okp, W2, r2 = check(R['prime_sep_cols'], X, mod, t)
    rat.append(W1/R['sc'])
    print(X, tt, 'prime pairs', ok1, 'W=%.2f r=%d (claim %.2f)' % (W1, r1, R['prime']), '| +sep pairs', ok2, 'placement', okp,
          'W=%.2f r=%d (claim %.2f)' % (W2, r2, R['prime_sep']), 'W_prime/W_sc_ub=%.3f' % (W1/R['sc']))
print('W_prime/W_sc(ub) range', min(rat), max(rat))
