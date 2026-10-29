"""Theory-lane spot checks (own code, mod Q1 and Q2 unless stated):
(a) the round-27 random boxes (n=6 seeds 100-114, n=7 seeds 1-12): h(sigma), rigidity (R), reach_pair; theory claims
    (R) holds on exactly 3 of the 27 and then reach = 2h.
(b) linkage identity H_A(e) = T_d(e) - H_B(2d-e), B = cube minus A, on syl8 L=2, L=3, sdc7 L=2.
(c) M8 with the 14 points |sigma| = 3 removed: reach_pair (theory: 7 = zeta_pair).
"""
import sys
import numpy as np
import flint
sys.path.insert(0, '/tmp/claude-0/-home-user-Jarnik/ba4929a4-8665-5790-985a-3f09130a14b9/scratchpad/round28/verify_checks')
from vlib import *

def Hof(S, d, p):
    by = std_monos(d, 2 * d)
    b = least_blocks(S, by, p, 2 * d)
    H = [x.shape[0] for x in b]
    while H and H[-1] == 0: H.pop()
    return H

def reach_of(A, d, p):
    by = std_monos(d, 2 * d)
    H = Hof(A, d, p); rho = len(H) - 1
    idx = {tuple(k.tolist()): i for i, k in enumerate(A)}
    sh = [least_blocks(shifted(A, np.array(l), idx), by, p, rho) for l in lpair_classes(d)]
    for e in range(rho + 1):
        rows = [r for b in sh for r in b[e].tolist()]
        r = flint.nmod_mat(rows, p).rank() if rows else 0
        if r < H[e]: return e, (r, H[e])
    return None, None

def rigid(A):
    d = A.shape[1]; idx = set(tuple(k.tolist()) for k in A)
    sig = A.astype(int).sum(axis=1); h = int(sig.max())
    for k in A[sig == h]:
        for x in range(d):
            for y in range(d):
                if x != y:
                    t = k.astype(int).copy(); t[x] += 1; t[y] -= 1
                    if tuple(t.tolist()) in idx: return h, False
    return h, True

print('(a) round-27 random boxes', flush=True)
cnt = 0
for n, seeds in ((6, range(100, 115)), (7, range(1, 13))):
    for sd in seeds:
        pats, L = rand_profile(n, sd); A = box(n, pats, L); d = n - 1
        h, R = rigid(A)
        r1 = reach_of(A, d, Q1); r2 = reach_of(A, d, Q2)
        cnt += R
        print('  n', n, 'seed', sd, 'L', L, '|A|', len(A), 'h', h, '(R)', R, 'reach', r1[0], r2[0], 'fail', r1[1],
              ('reach == 2h' if r1[0] == 2 * h else 'reach != 2h') if R else '', flush=True)
print('  (R) holds on', cnt, 'of 27', flush=True)

print('(b) linkage identity', flush=True)
for spec in ('syl:3:2', 'syl:3:3', 'sdc:7:2'):
    n, pats, L, A = spec_box(spec); d = n - 1
    C = cube(d); Aset = set(tuple(k.tolist()) for k in A)
    B = np.array([k for k in C if tuple(k.tolist()) not in Aset], dtype=np.int8)
    T = [1]
    for _ in range(d): T = list(np.convolve(T, [1, 1, 1]))
    for p in (Q1,):
        HA = Hof(A, d, p); HB = Hof(B, d, p)
        HA += [0] * (2 * d + 1 - len(HA)); HB += [0] * (2 * d + 1 - len(HB))
        ok = all(HA[e] == T[e] - HB[2 * d - e] for e in range(2 * d + 1))
        print('  ', spec, '|B|', len(B), 'identity holds for all e:', ok, flush=True)

print('(c) M8 minus |sigma|=3', flush=True)
n, pats, L, A = spec_box('syl:3:2')
sig = A.astype(int).sum(axis=1)
A2 = A[np.abs(sig) <= 2]
for p in (Q1, Q2):
    print('  p', p, '|A\'|', len(A2), 'H', Hof(A2, 7, p), 'reach_pair', reach_of(A2, 7, p), flush=True)
