#!/usr/bin/env python3
"""Exact checks for Section 2 of hardness.md (the L^2 / additive-energy side).

(A) Difference multiplicity: for I = [X, X+L] with L <= c*sqrt(X), every nonzero d has
    #{(x,x') in I^2 : x^2 - x'^2 = d} <= 1 + floor(c^2).
(B) Two-interval energy: E(A,B) <= |A||B| + (1+floor(c^2)) min(|A|,|B|)^2 <= (2+c^2)|A||B|.
(C) One-interval sums are trivially bounded: #{(x,y) in I^2 : x^2+y^2 = n} <= 2(c^2/sqrt2 + 2).
(D) Rudin in the far regime: if a >= (qK/2)^(4/3) the AP {a+qn : 0<=n<K} has <= sqrt(6K) squares.
All arithmetic is exact integer arithmetic (numpy int64 with explicit overflow guards, or Python ints).
"""
import math, random, sys
import numpy as np

random.seed(20261009)
out = []
def log(s):
    print(s); out.append(s)

def diff_mult_max(X, L):
    xs = np.arange(X, X + L + 1, dtype=np.int64)
    sq = xs * xs
    assert sq[-1] < 2**62
    d = (sq[:, None] - sq[None, :]).ravel()
    d = d[d > 0]
    if d.size == 0:
        return 0
    _, cnt = np.unique(d, return_counts=True)
    return int(cnt.max())

def sum_mult_max(I, J):
    xs = np.arange(I[0], I[1] + 1, dtype=np.int64)
    ys = np.arange(J[0], J[1] + 1, dtype=np.int64)
    s = (xs[:, None] ** 2 + ys[None, :] ** 2).ravel()
    _, cnt = np.unique(s, return_counts=True)
    return int(cnt.max())

# (A)
log("(A) difference multiplicity on I=[X,X+floor(c sqrt X)]")
worst = 0.0
for trial in range(60):
    c = random.choice([0.5, 1.0, 1.5, 2.0, 3.0])
    X = random.randint(10**4, 4 * 10**6)
    L = int(c * math.isqrt(X))
    while L * L > c * c * X:
        L -= 1
    m = diff_mult_max(X, L)
    bound = 1 + math.floor(c * c)
    assert m <= bound, (X, L, c, m)
    worst = max(worst, m / bound)
    if trial < 6:
        log(f"  X={X} L={L} c={c}: max_d r(d) = {m} <= {bound}")
log(f"  60 random intervals: all within bound; worst ratio {worst:.3f}")

# (B) energy for random subsets of two different intervals
log("(B) two-interval energy E(A,B) <= |A||B| + (1+floor c^2) min(|A|,|B|)^2")
for trial in range(40):
    c = random.choice([1.0, 2.0])
    X = random.randint(10**5, 10**6)
    Y = random.randint(10**4, X)
    L = int(c * math.isqrt(Y))
    A = sorted(random.sample(range(X, X + L + 1), random.randint(5, L)))
    B = sorted(random.sample(range(Y, Y + L + 1), random.randint(5, L)))
    a2 = np.array(A, dtype=np.int64) ** 2
    b2 = np.array(B, dtype=np.int64) ** 2
    s = (a2[:, None] + b2[None, :]).ravel()
    _, cnt = np.unique(s, return_counts=True)
    E = int((cnt.astype(np.int64) ** 2).sum())
    bnd = len(A) * len(B) + (1 + math.floor(c * c)) * min(len(A), len(B)) ** 2
    assert E <= bnd <= (2 + c * c) * len(A) * len(B)
    if trial < 4:
        log(f"  |A|={len(A)} |B|={len(B)} c={c}: E={E} <= {bnd}")
log("  40 random pairs of sets: all within bound")

# (C) one-interval sums
log("(C) one-interval sum multiplicity <= 2(c^2/sqrt2+2)")
for trial in range(30):
    c = random.choice([1.0, 2.0, 3.0])
    X = random.randint(10**5, 2 * 10**6)
    L = int(c * math.isqrt(X))
    m = sum_mult_max((X, X + L), (X, X + L))
    bound = 2 * (c * c / math.sqrt(2) + 2)
    assert m <= bound, (X, L, m)
    if trial < 4:
        log(f"  X={X} L={L}: max r(n) = {m} <= {bound:.2f}")
log("  30 random intervals: all within bound")
# contrast: two far-apart intervals at a generic direction (this is the uniform theorem itself)
log("  contrast, I and J far apart (generic direction): max r(n) over random boxes")
best = 0
for trial in range(30):
    X = random.randint(10**5, 10**6)
    Y = random.randint(X // 3, (2 * X) // 3)
    L = int(1.0 * math.isqrt(int(math.isqrt(X * X + Y * Y))))
    m = sum_mult_max((X, X + L), (Y, Y + L))
    best = max(best, m)
log(f"  30 random generic boxes of side sqrt(R): max r(n) found = {best}")

# (D) Rudin far regime
log("(D) squares in APs with a >= (qK/2)^(4/3): count <= sqrt(6K)")
worst = 0.0
for trial in range(3000):
    K = random.randint(2, 400)
    q = random.randint(1, 10**4)
    amin = math.ceil((q * K / 2) ** (4 / 3)) + 1
    # pick a so that the AP contains at least one square: a = y^2 - q*n0
    y = random.randint(math.isqrt(amin) + 2, math.isqrt(amin) * 3 + 5)
    n0 = random.randint(0, K - 1)
    a = y * y - q * n0
    if a < (q * K / 2) ** (4 / 3):
        continue
    lo = math.isqrt(a - 1) + 1
    hi = math.isqrt(a + q * (K - 1))
    cnt = sum(1 for t in range(lo, hi + 1) if (t * t - a) % q == 0)
    assert cnt <= math.sqrt(6 * K), (a, q, K, cnt)
    worst = max(worst, cnt / math.sqrt(6 * K))
log(f"  random APs: all within bound; worst count/sqrt(6K) = {worst:.3f}")

with open(__file__.replace('.py', '_output.txt'), 'w') as f:
    f.write('\n'.join(out) + '\n')
