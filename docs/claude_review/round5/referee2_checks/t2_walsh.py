r"""Referee-2 check T2: Remark 2.3 of construct.md (Walsh cores of order >= 16 have no skew pairing).

For H(x,a) = (-1)^(x.a) the skew-pairing minor condition for x, y != x_inf reads
    (-1)^(x.pi(x) + y.pi(y)) != (-1)^(x.pi(y) + y.pi(x))   <=>   (x+y).(pi(x)+pi(y)) = 1 (mod 2).
(1) We verify this equivalence directly on the matrix for t = 3, 4 (all x, y, all label pairs).
(2) Relaxed problem used by the restriction argument: injective p : F_2^t \ {0} -> F_2^t (0 allowed) with
    (x+y).(p(x)+p(y)) = 1 for all distinct x, y.  We COUNT all solutions by a bitset DFS with forward checking
    (a different implementation and variable order from construct's and the previous referee's scripts).
(3) Restriction lemma, checked by brute force on t = 5 sample data: for x, y in a subspace U, the quantity
    (x+y).(pi(x)+pi(y)) depends only on pi(x)|_U, pi(y)|_U.
Also the translation reduction: H(x + x0, a) = (-1)^(x0.a) H(x, a) (a column sign change), so x_inf = 0 w.l.o.g.
"""
import itertools, random, time
import numpy as np

def dot(a, b):
    return bin(a & b).count("1") & 1

t0 = time.time()
# (1) equivalence of the minor condition
for t in (3, 4):
    M = 1 << t
    H = np.array([[(-1) ** dot(x, a) for a in range(M)] for x in range(M)])
    for x, y in itertools.combinations(range(M), 2):
        for a in range(M):
            for b in range(M):
                lhs = H[x, a] * H[y, b] != H[x, b] * H[y, a]
                rhs = dot(x ^ y, a ^ b) == 1
                assert lhs == rhs
print("(1) minor condition <=> (x+y).(a+b) = 1 verified for t = 3, 4 (all x, y, a, b)")

# (2) count relaxed solutions
def count(t, allow_zero, bijective_nonzero=False):
    M = 1 << t
    xs = list(range(1, M))
    # order variables by popcount descending (different from natural order)
    xs.sort(key=lambda x: (-bin(x).count("1"), x))
    vals = list(range(0 if allow_zero else 1, M))
    nsol = [0]; nodes = [0]
    assigned = []
    def dfs(i, used):
        nodes[0] += 1
        if i == len(xs):
            nsol[0] += 1
            return
        x = xs[i]
        for v in vals:
            if used >> v & 1:
                continue
            ok = True
            for (y, w) in assigned:
                if dot(x ^ y, v ^ w) != 1:
                    ok = False; break
            if ok:
                assigned.append((x, v))
                dfs(i + 1, used | (1 << v))
                assigned.pop()
    dfs(0, 0)
    return nsol[0], nodes[0]

for t in (2, 3):
    s, nd = count(t, allow_zero=False)
    print("(2) t=%d: bijections F_2^t\\0 -> F_2^t\\0 with the condition (core skew pairings, x_inf = 0): %d  (nodes %d)" % (t, s, nd))
    s, nd = count(t, allow_zero=True)
    print("    t=%d: injections F_2^t\\0 -> F_2^t (0 allowed): %d  (nodes %d)" % (t, s, nd))
s, nd = count(4, allow_zero=True)
print("(2) t=4: injections F_2^4\\0 -> F_2^4 (0 allowed) with the condition: %d solutions (nodes %d)" % (s, nd))
assert s == 0

# (3) restriction lemma, brute-force sanity on random data for t = 5, U = span of 4 random independent vectors
rng = random.Random(5)
t = 5; M = 1 << t
for trial in range(200):
    while True:
        basis = [rng.randrange(1, M) for _ in range(4)]
        span = {0}
        for bvec in basis:
            span |= {s ^ bvec for s in span}
        if len(span) == 16:
            break
    U = sorted(span)
    pi = {x: rng.randrange(M) for x in U}
    # restriction of functional pi(x) to U: the vector of values u -> u.pi(x) on U
    def restr(v):
        return tuple(dot(u, v) for u in U)
    for x, y in itertools.combinations([u for u in U if u], 2):
        val = dot(x ^ y, pi[x] ^ pi[y])
        # recompute from restrictions only: (x+y) is in U, so value = restr(pi x)[idx(x+y)] + restr(pi y)[idx(x+y)]
        idx = U.index(x ^ y)
        val2 = (restr(pi[x])[idx] + restr(pi[y])[idx]) % 2
        assert val == val2
print("(3) restriction lemma: (x+y).(pi x + pi y) depends only on restrictions to U (200 random t=5 trials): OK")
# translation reduction
for t in (3, 4):
    M = 1 << t
    for x0 in range(M):
        for x in range(M):
            for a in range(M):
                assert (-1) ** dot(x ^ x0, a) == (-1) ** dot(x0, a) * (-1) ** dot(x, a)
print("translation reduction verified for t = 3, 4")
print("T2 DONE (%.1fs)" % (time.time() - t0))
