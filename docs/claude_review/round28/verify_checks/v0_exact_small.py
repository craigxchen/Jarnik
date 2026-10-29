"""Validate the mod-p engine against an exact-over-Q computation (fmpq rref), all degrees, BOTH signs of every
L_pair shift, on small slab boxes.  Also checks M_e(A_{-l}) = M_e(A_l) exactly."""
import sys, random
import numpy as np
import flint
sys.path.insert(0, '/tmp/claude-0/-home-user-Jarnik/ba4929a4-8665-5790-985a-3f09130a14b9/scratchpad/round28/verify_checks')
from vlib import *

def blocks_Q(S, by, emax):
    monos = [a for e in range(emax + 1) for a in by[e]]
    st = [0]
    for e in range(emax + 1): st.append(st[-1] + len(by[e]))
    E = evalmat(S, monos).astype(int).tolist()
    R, r = flint.fmpq_mat(E).rref()
    piv = []; j = 0
    for i in range(r):
        while R[i, j] == 0: j += 1
        piv.append(j); j += 1
    out = []
    for e in range(emax + 1):
        rows = [[R[i, c] for c in range(st[e], st[e + 1])] for i in range(r) if st[e] <= piv[i] < st[e + 1]]
        out.append(rows)
    return out

def rkQ(rows):
    return flint.fmpq_mat(rows).rank() if rows else 0

rng = random.Random(2828)
boxes = [('syl:2:1',), ('syl:3:1',)]
specs = []
for t in range(14):
    n = rng.choice([4, 5])
    r = min(rng.randint(n - 1, 2 * n), 2 ** (n - 1) - 1)
    pats = set()
    while len(pats) < r:
        S = [x for x in range(n) if rng.random() < 0.5]
        if 0 < len(S) < n: pats.add(canon(n, S))
    specs.append((n, sorted(pats), rng.choice([1, 2])))
specs.append((4, [], 1))  # full cube d=3
n2, p2 = syl(3); specs.append((n2, p2, 1))
bad = 0
for (n, pats, L) in specs:
    A = box(n, pats, L); d = n - 1
    by = std_monos(d, 2 * d)
    bQ = blocks_Q(A, by, 2 * d)
    HQ = [len(b) for b in bQ]
    while HQ[-1] == 0: HQ.pop()
    rho = len(HQ) - 1
    idx = {tuple(k.tolist()): i for i, k in enumerate(A)}
    rows_out = []
    for p in (Q1, Q2):
        bp = least_blocks(A, by, p, 2 * d)
        Hp = [b.shape[0] for b in bp][:rho + 1]
        if Hp != HQ: bad += 1; print('H mismatch', p, Hp, HQ)
    shiftsQ = {}
    for l in lpair_classes(d):
        for sgn in (1, -1):
            ll = np.array(l) * sgn
            Al = shifted(A, ll, idx)
            shiftsQ[(l, sgn)] = blocks_Q(Al, by, rho) if len(Al) else [[] for _ in range(rho + 1)]
    reachQ = None; dimsQ = []
    for e in range(rho + 1):
        allrows = [r for k in shiftsQ for r in shiftsQ[k][e]]
        plus = [r for (l, s) in shiftsQ if s == 1 for r in shiftsQ[(l, s)][e]]
        dq = rkQ(allrows); dqp = rkQ(plus)
        dimsQ.append((dq, HQ[e]))
        if dq != dqp: bad += 1; print('sign mismatch', e, dq, dqp)
        # sign check per shift
        for l in lpair_classes(d):
            a = shiftsQ[(l, 1)][e]; b = shiftsQ[(l, -1)][e]
            if rkQ(a) != rkQ(b) or rkQ(a + b) != rkQ(a): bad += 1; print('M_e(A_l) != M_e(A_-l)', l, e)
        if reachQ is None and dq < HQ[e]: reachQ = e
    # mod p, own driver logic inline (one sign)
    for p in (Q1, Q2):
        reachp = None; dimsp = []
        for e in range(rho + 1):
            rows = []
            for l in lpair_classes(d):
                Al = shifted(A, np.array(l), idx)
                if len(Al) > e: rows += least_blocks(Al, by, p, e)[e].tolist()
            dp = flint.nmod_mat(rows, p).rank() if rows else 0
            dimsp.append((dp, HQ[e]))
            if reachp is None and dp < HQ[e]: reachp = e
        if dimsp != dimsQ: bad += 1; print('dims mismatch', p, dimsp, dimsQ)
    print('n', n, 'cuts', len(pats), 'L', L, '|A|', len(A), 'H', HQ, 'reachQ', reachQ, 'dimsQ', dimsQ, flush=True)
print('disagreements', bad)
