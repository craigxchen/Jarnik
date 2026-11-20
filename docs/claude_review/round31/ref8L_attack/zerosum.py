"""Zero-sum subset counts for vectors z_j = u_j * (1, n_j, ..., n_j^{r-1}) (moment curve), or for
arbitrary pairwise non-parallel vectors when r = 2 (the moment-curve condition is then vacuous,
node infinity <-> direction (0,1)).

Dictionary (proved in the report): k nodes, r = k - d - 1.  The polynomials of degree <= d taking
two values at every node correspond to 0/1 vectors x with sum_j x_j z_j = 0, z_j = u_j m(n_j),
for suitable weights u_j; B = number of such x (x = 0 included), exponent E = k/d.
"""
import random, sys, itertools
from fractions import Fraction as Fr

def count_zero_sums(Z):
    k = len(Z); r = len(Z[0])
    # meet in the middle on exact integer vectors
    h = k // 2
    left = {}
    for m in range(1 << h):
        s = tuple(sum(Z[j][c] for j in range(h) if (m >> j) & 1) for c in range(r))
        left.setdefault(s, []).append(m)
    cnt = 0; sols = []
    for m in range(1 << (k - h)):
        s = tuple(-sum(Z[h + j][c] for j in range(k - h) if (m >> j) & 1) for c in range(r))
        if s in left:
            for a in left[s]:
                cnt += 1; sols.append(a | (m << h))
    return cnt, sols

def general_position(Z):
    import sympy
    r = len(Z[0])
    for S in itertools.combinations(range(len(Z)), r):
        if sympy.Matrix([Z[j] for j in S]).det() == 0: return False
    return True

def moment_vec(u, n, r):
    if n == 'inf': return [0] * (r - 1) + [u]
    return [u * n ** p for p in range(r)]

def to_int(Z):
    from math import lcm
    den = 1
    for z in Z:
        for x in z: den = lcm(den, Fr(x).denominator)
    return [[int(Fr(x) * den) for x in z] for z in Z]

if __name__ == "__main__":
    r = int(sys.argv[1]); k = int(sys.argv[2]); trials = int(sys.argv[3]); seed = int(sys.argv[4])
    random.seed(seed)
    npool = list(range(-6, 7)) + ['inf']
    upool = [x for x in range(-6, 7) if x]
    best = (0,)
    for t in range(trials):
        nodes = random.sample(npool, k)
        # weights: choose to create zero sums -- random small integers
        us = [random.choice(upool) for _ in range(k)]
        Z = to_int([moment_vec(Fr(u), n if n == 'inf' else Fr(n), r) for u, n in zip(us, nodes)])
        c, sols = count_zero_sums(Z)
        if c > best[0]:
            best = (c, nodes, us)
    print(f"r={r} k={k} d={k-r-1} E={k/(k-r-1):.4f} best B={best[0]} nodes={best[1]} u={best[2]}")
