"""Referee: recompute the exact independent mean E D and aggregate §2.1 statistics from genuine_trials.jsonl."""
import json, math
import numpy as np
from sympy import primerange
F = '/tmp/claude-0/-home-user-Jarnik/ba4929a4-8665-5790-985a-3f09130a14b9/scratchpad/round7/p3data_checks/genuine_trials.jsonl'
def mQ(q, a): return (q - 1 if q % 4 == 1 else q + 1) * q**(a-1)
def ED(X):
    s = 0
    for q in primerange(3, X+1):
        a = 1
        while q**a <= X:
            s += 2*math.log(q)/mQ(q, a); a += 1
    return s
recs = [json.loads(l) for l in open(F)]
for M in [12, 20, 24, 32, 44, 60, 72, 84]:
    for k, Xf in (('u1_X2', 2), ('u1_X4', 4)):
        R = [r for r in recs if r['M'] == M and k in r]
        L = math.log(M); H = L*L/math.log(L)
        mean = np.mean([r[k]['mean'] for r in R]); dmax = np.median([r[k]['max'] for r in R])
        lam = np.mean([r[k]['lam'] for r in R]); mrow = np.mean([r[k]['maxrow'] for r in R])
        pen = np.median([(M-1)*r[k]['max']/r[k]['lam'] for r in R]); pen2 = (M-1)*dmax/lam
        print(M, Xf*M, len(R), 'ED %.3f mean %.3f Dmax med %.1f Dmax/L %.2f Dmax/H %.2f lam/ML %.3f maxrow/ML %.3f pen(med of ratio) %.2f pen(ratio of agg) %.2f' % (
            ED(Xf*M), mean, dmax, dmax/L, dmax/H, lam/(M*L), mrow/(M*L), pen, pen2))
