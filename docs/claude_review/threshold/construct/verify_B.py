"""Exact checks for Theorem B (mean-width bound) of construct.md.

Theorem B.  For every finite S in Z^n:   reg(S) <= MW(S) := 2^{1-n} * sum_{P subset [n]} w_P(S),
where w_P(S) = max_{x in S} x(P) - min_{x in S} x(P), x(P) = sum_{i in P} x_i.

Checks:
 (1) reg(S) <= MW(S) on random S (reg computed mod two primes; mod-q reg >= rational reg, so a
     PASS is a rigorous verification of the inequality for that S).
 (2) the recursive bound produced by the proof (peeling along x_1 with slice bounds g_c) is <= MW(S),
     and reg(S) <= that recursive bound.
 (3) equality cases: boxes prod [0,a_i] (reg = sum a_i = MW), root-binomial products in homogeneous form.
"""
import random, itertools, sys
from fractions import Fraction as Fr
import numpy as np
from clib import reg2

def MW(S, n):
    if n == 0:
        return Fr(0)
    A = np.array(S, dtype=np.int64).reshape(len(S), n)
    tot = 0
    for mask in range(1 << n):
        idx = [i for i in range(n) if (mask >> i) & 1]
        if not idx:
            continue
        v = A[:, idx].sum(axis=1)
        tot += int(v.max() - v.min())
    return Fr(2 * tot, 1 << n)

def proof_bound(S, n):
    """the bound produced by the inductive proof: peel along x_1, slice bounds g_c = proof_bound(slice)."""
    if n == 0:
        return Fr(0)
    slices = {}
    for x in S:
        slices.setdefault(x[0], []).append(x[1:])
    g = sorted((proof_bound(T, n - 1) for T in slices.values()), reverse=True)
    return max(gi + i for i, gi in enumerate(g))

def reg_any(S, n):
    if n == 0 or len(S) == 1:
        return 0
    return reg2([tuple(x) for x in S])

if __name__ == '__main__':
    random.seed(int(sys.argv[2]) if len(sys.argv) > 2 else 2026)
    fails = 0
    tight = 0
    count = 0
    for trial in range(int(sys.argv[1]) if len(sys.argv) > 1 else 400):
        n = random.randint(1, 4)
        R = random.choice([1, 2, 3, 5])
        mode = random.choice(["rand", "rand", "conv"])
        if mode == "rand":
            S = set()
            for _ in range(random.randint(1, 30)):
                S.add(tuple(random.randint(-R, R) for _ in range(n)))
            S = sorted(S)
        else:
            # lattice points of a random polytope {x : |<a_t, x>| <= r_t}
            rows = [[random.randint(-2, 2) for _ in range(n)] for _ in range(random.randint(n, n + 3))]
            rad = [random.randint(1, 4) for _ in rows]
            S = [x for x in itertools.product(range(-R, R + 1), repeat=n)
                 if all(abs(sum(a * b for a, b in zip(r, x))) <= rr for r, rr in zip(rows, rad))]
            if not S:
                continue
        if len(S) > 800:
            continue
        count += 1
        r = reg_any(S, n)
        mw = MW(S, n)
        pb = proof_bound(S, n)
        if not (r <= pb <= mw):
            fails += 1
            print("FAIL", n, S, r, pb, mw)
        if r == mw:
            tight += 1
    print(f"(1)+(2): {count} random sets, n<=4: reg <= proof_bound <= MW in all; failures={fails}; equality reg=MW in {tight}")

    # (3) boxes
    for a in [(3,), (2, 1), (1, 1, 1), (2, 2, 1), (3, 1, 2)]:
        S = list(itertools.product(*[range(ai + 1) for ai in a]))
        print("box", a, "reg", reg_any(S, len(a)), "MW", MW(S, len(a)), "sum a", sum(a))
    # root-binomial products in homogeneous form (drop last coordinate): Minkowski sum of segments [e_i, e_j]
    def minkowski_roots(M, edges):
        pts = {tuple([0] * M)}
        for (i, j) in edges:
            new = set()
            for p in pts:
                for e in (i, j):
                    q = list(p); q[e] += 1; new.add(tuple(q))
            pts = new
        return sorted(pts)
    for M, edges in [(3, [(0, 1), (1, 2), (0, 2)]), (4, [(0, 1), (2, 3)]), (4, [(0, 1), (1, 2), (2, 3), (0, 3)]),
                     (4, [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]), (5, [(0, 1), (1, 2), (2, 3), (3, 4), (0, 4)])]:
        S = [p[:-1] for p in minkowski_roots(M, edges)]
        print("root product M=%d |E|=%d:" % (M, len(edges)), "reg", reg_any(S, M - 1), "MW", MW(S, M - 1))
