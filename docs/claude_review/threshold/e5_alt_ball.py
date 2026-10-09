# E5: alternating (S_M-antisymmetric) measures on the BALL A_{s,D} versus on the shell S_{s,D}.
# An antisymmetric measure mu = sum_rep c(rep) Alt(delta_rep) has order C(M,2) + e, where e is the
# largest integer such that some nonzero c''' on the sorted reps (distinct entries) annihilates all
# symmetric polynomials of degree < e, i.e. polynomials in p_2..p_M (p_1 = s fixed) of weighted
# degree < e.  Computed mod p (rank mod p <= rank over Q, so e_p >= e_Q: an UPPER bound on the order).
import sys, numpy as np
from math import comb
from taulib import P1, P2, pivot_columns_mod

def reps(M, s, D, shell_only):
    out = []
    def rec(prefix, maxv, rem_sum, rem_l1, left):
        if left == 0:
            if rem_sum == 0 and (rem_l1 == 0 or (not shell_only and rem_l1 % 2 == 0)):
                out.append(tuple(prefix))
            return
        for v in range(min(maxv, rem_l1), -rem_l1 - 1, -1):
            if abs(v) > rem_l1: continue
            rec(prefix + [v], v - 1, rem_sum - v, rem_l1 - abs(v), left - 1)
    rec([], D, s, D, M)
    return out

def wmonomials(M, e):
    parts = list(range(2, M + 1)); out = []
    def rec(i, rem, cur):
        if i == len(parts): out.append(tuple(cur)); return
        w = parts[i]
        for a in range(rem // w + 1): rec(i + 1, rem - a * w, cur + [a])
    rec(0, e, []); return out

def e_max(R, M, p):
    n = len(R)
    ps = [[pow(sum(x ** j for x in r), 1, p) % p for j in range(2, M + 1)] for r in R]
    e = 0
    while True:
        mons = wmonomials(M, e)
        A = np.zeros((n, len(mons)), dtype=np.int64)
        for i, pv in enumerate(ps):
            for jj, a in enumerate(mons):
                v = 1
                for b, ex in zip(pv, a):
                    if ex: v = v * pow(b, ex, p) % p
                A[i, jj] = v
        if len(pivot_columns_mod(A, p)) == n: return e
        e += 1

M = int(sys.argv[1]); extra = int(sys.argv[2])
v = list(range(-(M // 2), M - M // 2)); tmin = sum(abs(x) for x in v)
best = {'shell': (0,), 'ball': (0,)}
for D in range(tmin, tmin + extra + 1):
    for s in range(0, 4):
        if (D - s) % 2: continue
        for mode in ('shell', 'ball'):
            R = reps(M, s, D, mode == 'shell')
            if not R: continue
            e1 = e_max(R, M, P1); e2 = e_max(R, M, P2)
            order = comb(M, 2) + e1
            flag = '' if e1 == e2 else ' PRIME-MISMATCH'
            if order / D > best[mode][0]: best[mode] = (order / D, D, s, len(R), order)
            print(f"M={M} {mode:5s} D={D} s={s} reps={len(R)} e={e1} order={order} ratio={order/D:.4f}{flag}", flush=True)
print("BEST", M, best, flush=True)
