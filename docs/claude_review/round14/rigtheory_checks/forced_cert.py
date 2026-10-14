"""Forced sets for K=2..7, certified complete: scan every prime q <= 4A_K+5 (Theorem 1(b) of round10 L4K.md
with a greedy B_h set, verified here), above which every q is free.  Usage: python3 forced_cert.py p1 p2 ..."""
import sys, itertools
from rlib14 import *
from forced import forced_set
def greedy_Bh(r, h):
    S = []
    x = 0
    while len(S) < r:
        x += 1
        T = S + [x]
        sums = set(); ok = True
        for comb in itertools.combinations_with_replacement(range(len(T)), h):
            s = sum(T[i] for i in comb)
            if s in sums: ok = False; break
            sums.add(s)
        if ok: S.append(x)
    return S
ps = [int(x) for x in sys.argv[1:]]; r = len(ps)
res = {}
for K in range(2, 8):
    h = K // 2
    beta = list(range(1, r + 1)) if h <= 1 else greedy_Bh(r, h)
    B = max(beta)
    AK = K * (2 * (K + 1) ** 2 * B + 2 * K + 1)
    qT = 4 * AK + 5
    F = forced_set(ps, K, qT)
    res[K] = F
    print(f"K={K} h={h} B_h set={beta} A_K={AK} scanned q<={qT}: #forced={len(F)} q*={max(F)}", flush=True)
    print("  forced:", F, flush=True)
