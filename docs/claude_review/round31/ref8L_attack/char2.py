"""Characteristic-2 check: over F_q[T], q = 2^m, nodes = all of F_q, unit weights.
For each subset A of F_q put Y_A(T) = sum_{n in A} prod_{l != n} (T - l).  Then Y_A(n) = [n in A]
(since prod_{l != n}(n - l) = -1 = 1), so the Y_A take two values at every node, and
deg Y_A <= q - 1 - r where r = number of vanishing moments sum_{n in A} n^s, s = 0..r-1.
Claim: the 2q sets {empty, F_q, affine F_2-hyperplanes} have r = q/2 - 1, so deg Y_A <= q/2
(exponent q/(q/2) = 2) with B = 2q points.  We list ALL subsets with deg <= q/2 for q = 8, 16."""
import itertools, sys
MOD = {3: 0b1011, 4: 0b10011}
def mul(a, b, m):
    r = 0
    while b:
        if b & 1: r ^= a
        b >>= 1; a <<= 1
        if a >> m: a ^= MOD[m]
    return r
def pmul(P, Q, m):
    R = [0] * (len(P) + len(Q) - 1)
    for i, a in enumerate(P):
        for j, b in enumerate(Q): R[i + j] ^= mul(a, b, m)
    return R
def deg(P):
    for i in range(len(P) - 1, -1, -1):
        if P[i]: return i
    return -1
for m in (3, 4):
    q = 2 ** m
    basis = []
    for n in range(q):
        P = [1]
        for l in range(q):
            if l != n: P = pmul(P, [l, 1], m)   # (T + l) = (T - l) in char 2
        basis.append(P)
    good = []
    for mask in range(1 << q):
        Y = [0] * q
        for n in range(q):
            if mask >> n & 1:
                Y = [a ^ b for a, b in zip(Y, basis[n])]
        if deg(Y) <= q // 2: good.append(mask)
    # count affine hyperplanes
    hyper = set()
    for f in range(1, q):        # linear functional n -> parity(f & n) (F_2-linear)
        for c in (0, 1):
            hyper.add(sum(1 << n for n in range(q) if bin(f & n).count('1') % 2 == c))
    expected = hyper | {0, (1 << q) - 1}
    print(f"q={q}: subsets with deg Y_A <= q/2: {len(good)}; equals {{0, F_q, affine hyperplanes}}: {set(good) == expected} (size {len(expected)})")
