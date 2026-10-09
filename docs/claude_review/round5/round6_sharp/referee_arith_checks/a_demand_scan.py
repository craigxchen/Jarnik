"""Pair-demand scan over many M (evidence; complements a_assign.py Part D).
D(d) = sum_{q < 4M, q <= X, p'(q) | d} log q * min(a*(q), A_q, 1 + v_q(d)); X = 3M (A = 1) and X = 3M^2.
Reports max over M in ranges of max_d D(d)/log M."""
import math, numpy as np
LIM = 2 ** 22
s = np.ones(LIM + 1, dtype=bool); s[:2] = False
for i in range(2, int(LIM ** 0.5) + 1):
    if s[i]:
        s[i * i::i] = False
pr = np.nonzero(s)[0]
pp = {3: 2, 5: 2, 7: 3}
for k in range(3, 21):
    B = pr[(pr >= 2 ** k) & (pr < 2 ** (k + 1)) & (pr > 7)]
    T = pr[(pr >= 2 ** (k - 2)) & (pr < 2 ** (k - 1))]
    n, t = len(B), len(T)
    for i, q in enumerate(B.tolist()):
        pp[q] = int(T[(i * t) // n])
def demand(M, X):
    Dv = np.zeros(M)
    for q in pr[(pr > 2) & (pr < 4 * M)].tolist():
        if q > X: break
        p1 = pp[q]
        if p1 >= M: continue
        lq = math.log(q); Aq = 1
        while Aq < 64 and q ** (Aq + 1) <= X: Aq += 1
        ast = 1
        while p1 * q ** (ast - 1) < M: ast += 1
        top = min(ast, Aq)
        Dv[p1::p1] += lq
        a = 2
        while a <= top:
            st = p1 * q ** (a - 1)
            if st >= M: break
            Dv[st::st] += lq; a += 1
    Dv[0] = 0
    return Dv.max()
for lo, hi, step in [(20, 200, 1), (200, 2000, 1), (2000, 20000, 97)]:
    best = (0, None)
    best2 = (0, None)
    for M in range(lo, hi, step):
        v = demand(M, 3 * M) / math.log(M)
        if v > best[0]: best = (v, M)
        v2 = demand(M, 3 * M * M) / math.log(M)
        if v2 > best2[0]: best2 = (v2, M)
    print(f"M in [{lo},{hi}) step {step}: max demand/log M = {best[0]:.3f} at M = {best[1]} (A=1);  {best2[0]:.3f} at M = {best2[1]} (A=2)")
