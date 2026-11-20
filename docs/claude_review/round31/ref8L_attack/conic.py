"""Moment-curve zero-sum search, two-block seeding.
Nodes n_1..n_k (rationals or 'inf'), r = k - d - 1.  Partition nodes into blocks; on each block of
size r+1 the weights are forced (up to scale) to the dual weights 1/prod(n_j - n_l); the relative
scales are chosen among the values that make some further subset zero-sum.  Counts all zero-sum
subsets exactly."""
import sys, itertools
from fractions import Fraction as Fr
from collections import Counter

def mvec(n, r):
    if n == 'inf': return tuple([Fr(0)] * (r - 1) + [Fr(1)])
    return tuple(Fr(n) ** p for p in range(r))

def dual_weights(block):
    # weights w_j with sum_j w_j n_j^p = 0, p < len(block)-1; infinity handled by projective formula
    fin = [n for n in block if n != 'inf']
    w = []
    for n in block:
        if n == 'inf':
            # coefficient of top power: for finite nodes F, sum_F n^p / prod_{l in F, l!=j}(n_j-n_l) = 0 for p<|F|-1
            # and = 1 for p = |F|-1, so the infinity weight is -1.
            w.append(Fr(-1))
        else:
            pr = Fr(1)
            for l in fin:
                if l != n: pr *= (Fr(n) - Fr(l))
            w.append(1 / pr)
    return w

def count(Z):
    k = len(Z); r = len(Z[0]); h = k // 2; left = Counter()
    for m in range(1 << h):
        left[tuple(sum((Z[j][c] for j in range(h) if m >> j & 1), Fr(0)) for c in range(r))] += 1
    c = 0
    for m in range(1 << (k - h)):
        s = tuple(-sum((Z[h + j][c2] for j in range(k - h) if m >> j & 1), Fr(0)) for c2 in range(r))
        c += left.get(s, 0)
    return c

def two_block(nodes, r):
    k = len(nodes); best = (0, None)
    for S in itertools.combinations(range(k), r + 1):
        T = [j for j in range(k) if j not in S]
        if len(T) != r + 1 or S[0] != 0: continue
        wS = dual_weights([nodes[j] for j in S]); wT = dual_weights([nodes[j] for j in T])
        A = {j: tuple(wS[i] * x for x in mvec(nodes[j], r)) for i, j in enumerate(S)}
        Bv = {j: tuple(wT[i] * x for x in mvec(nodes[j], r)) for i, j in enumerate(T)}
        ratios = set()
        for a in range(1, 1 << (r + 1)):
            va = [sum((A[S[i]][c] for i in range(r + 1) if a >> i & 1), Fr(0)) for c in range(r)]
            for b in range(1, 1 << (r + 1)):
                vb = [sum((Bv[T[i]][c] for i in range(r + 1) if b >> i & 1), Fr(0)) for c in range(r)]
                # lambda*va + vb = 0
                nz = [c for c in range(r) if va[c] != 0]
                if not nz: continue
                lam = -vb[nz[0]] / va[nz[0]]
                if lam != 0 and all(lam * va[c] + vb[c] == 0 for c in range(r)): ratios.add(lam)
        for lam in ratios:
            Z = [None] * k
            for j in S: Z[j] = tuple(lam * x for x in A[j])
            for j in T: Z[j] = Bv[j]
            c = count(Z)
            if c > best[0]: best = (c, (S, lam))
        if not ratios and best[0] < 4: best = (4, (S, None))
    return best

if __name__ == "__main__":
    r = int(sys.argv[1]); k = int(sys.argv[2]); lo, hi = int(sys.argv[3]), int(sys.argv[4])
    pool = [Fr(x) for x in range(lo, hi + 1)] + ['inf']
    best = (0,)
    for nodes in itertools.combinations(pool, k):
        c, info = two_block(list(nodes), r)
        if c > best[0]:
            best = (c, nodes, info); print(f"r={r} k={k} d={k-r-1} E={k/(k-r-1):.3f} B={c} nodes={[str(x) for x in nodes]} info={info}", flush=True)
