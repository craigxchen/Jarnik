"""Theory Thm 2.1 sanity check (own code): for products of small slab boxes, H convolves, reach(L_pair) adds,
zeta_pair adds.  Product profile: block-1 patterns (side omitting the base) on coords 0..d1-1, block-2 on d1..d1+d2-1,
shared base.  Reach mod Q1 and Q2."""
import sys, random
import numpy as np
import flint
sys.path.insert(0, '/tmp/claude-0/-home-user-Jarnik/ba4929a4-8665-5790-985a-3f09130a14b9/scratchpad/round28/verify_checks')
from vlib import *
from zono_v import functionals, candidates, lp_dual, DFS

def side(n, S):
    return tuple(sorted(set(range(n - 1)) - set(S))) if (n - 1) in S else tuple(S)

def reach_of(A, d, p):
    by = std_monos(d, 2 * d)
    bA = least_blocks(A, by, p, 2 * d)
    H = [b.shape[0] for b in bA]
    while H[-1] == 0: H.pop()
    rho = len(H) - 1
    idx = {tuple(k.tolist()): i for i, k in enumerate(A)}
    sh = [least_blocks(shifted(A, np.array(l), idx), by, p, rho) for l in lpair_classes(d)]
    for e in range(rho + 1):
        rows = [r for b in sh for r in b[e].tolist()]
        r = flint.nmod_mat(rows, p).rank() if rows else 0
        if r < H[e]: return H, e
    return H, None

def zpair(n, pats, A):
    d = n - 1
    F, W = functionals(n, pats, A)
    C = candidates(d, F, W, pair_only=True)
    Y, m, B, cv = lp_dual(F, W, C)
    z, w = DFS(A, C, F, Y, m).run()
    return z

rng = random.Random(77)
def rprof(n):
    r = min(rng.randint(n - 1, 2 * n), 2 ** (n - 1) - 1); pats = set()
    while len(pats) < r:
        S = [x for x in range(n) if rng.random() < 0.5]
        if 0 < len(S) < n: pats.add(canon(n, S))
    return sorted(pats)

tests = [((3, [(0, 1)], 1), (4, None, 1)), ((4, None, 1), (4, None, 2)), ((3, [(0, 1)], 1), (5, None, 2)),
         ((4, None, 2), (4, None, 1))]
for (n1, p1, L1), (n2, p2, L2) in tests:
    p1 = p1 if p1 is not None else rprof(n1)
    p2 = p2 if p2 is not None else rprof(n2)
    if L1 != L2:
        L2 = L1  # one common L (T=1 boxes have one L)
    d1, d2 = n1 - 1, n2 - 1; n = d1 + d2 + 1
    P = [side(n1, S) for S in p1] + [tuple(x + d1 for x in side(n2, S)) for S in p2]
    A1, A2 = box(n1, p1, L1), box(n2, p2, L2)
    A = box(n, [tuple(S) for S in P], L1)
    assert len(A) == len(A1) * len(A2)
    for p in (Q1, Q2):
        H1, r1 = reach_of(A1, d1, p); H2, r2 = reach_of(A2, d2, p); H, r = reach_of(A, d1 + d2, p)
        conv = list(np.convolve(H1, H2))
        print('p', p, 'H1', H1, 'H2', H2, 'H', H, 'conv ok', conv == H, 'reach', r1, r2, r, 'adds', r == r1 + r2, flush=True)
    z1, z2, z = zpair(n1, p1, A1), zpair(n2, p2, A2), zpair(n, P, A)
    print('  zeta_pair', z1, z2, z, 'adds', z == z1 + z2, flush=True)
