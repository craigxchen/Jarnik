"""Equality cases of Theorem B: find measures c with ord(c) = MW(supp c) on small random sets
(exact rational kernels) and test whether c is a product of root binomials (x_i - 1), (x_i - x_j)
times a monomial (dehomogenized coordinates; x_0 = 1 is the homogenizing coordinate)."""
import random, itertools, sys
from fractions import Fraction as Fr
from verify_B import MW

def monos(n, d):
    out = []
    for tot in range(d + 1):
        for c in itertools.combinations_with_replacement(range(n), tot):
            e = [0] * n
            for i in c: e[i] += 1
            out.append(tuple(e))
    return out

def kernel_order(S, n, m):
    """basis of {c : sum c_k k^a = 0 for |a| < m} over Q (exact)."""
    rows = []
    for a in monos(n, m - 1):
        rows.append([Fr(1) if True else 0 for _ in S])
        for j, k in enumerate(S):
            v = 1
            for i in range(n): v *= k[i] ** a[i]
            rows[-1][j] = Fr(v)
    # nullspace of rows (|rows| x |S|)
    A = [r[:] for r in rows]; ncol = len(S); piv = []; r = 0
    for c in range(ncol):
        p = next((i for i in range(r, len(A)) if A[i][c] != 0), None)
        if p is None: continue
        A[r], A[p] = A[p], A[r]
        inv = 1 / A[r][c]; A[r] = [x * inv for x in A[r]]
        for i in range(len(A)):
            if i != r and A[i][c] != 0:
                f = A[i][c]; A[i] = [x - f * y for x, y in zip(A[i], A[r])]
        piv.append(c); r += 1
        if r == len(A): break
    free = [c for c in range(ncol) if c not in piv]
    basis = []
    for fcol in free:
        v = [Fr(0)] * ncol; v[fcol] = Fr(1)
        for i, pc in enumerate(piv): v[pc] = -A[i][fcol]
        basis.append(v)
    return basis

def evalf(c, S, pt):
    tot = Fr(0)
    for ck, k in zip(c, S):
        if ck:
            t = ck
            for i in range(len(k)): t *= Fr(pt[i]) ** k[i]
            tot += t
    return tot

def strip_roots(c, S, n, rng):
    """repeatedly divide by root binomials; returns number of root factors found and whether the
    quotient is a monomial.  Divisibility by (x_a - x_b) (x_{-1}=1) tested at random rational points."""
    import copy
    # represent f as dict exponent->coef
    f = {k: ck for k, ck in zip(S, c) if ck}
    count = 0
    def vanishes_on(f, a, b):
        for _ in range(3):
            pt = [Fr(rng.randint(2, 50), rng.randint(1, 7)) for _ in range(n)]
            if b == -1: pt[a] = Fr(1)
            else: pt[a] = pt[b]
            val = sum(cf * prod_pow(pt, k) for k, cf in f.items())
            if val != 0: return False
        return True
    def prod_pow(pt, k):
        t = Fr(1)
        for i in range(n): t *= pt[i] ** k[i]
        return t
    def divide(f, a, b):
        # divide Laurent polynomial f by (x_a - x_b) (x_b := 1 if b == -1): synthetic division in x_a
        # write f = sum_j x_a^j g_j(other); quotient q with f = (x_a - y) q, y = x_b or 1
        g = {}
        for k, cf in f.items():
            g.setdefault(k[a], {})[k[:a] + (0,) + k[a+1:]] = cf
        js = sorted(g); q = {}; carry = {}
        # descending division
        rem = {j: dict(g[j]) for j in js}
        top = js[-1]; low = js[0]
        cur = {}
        for j in range(top, low, -1):
            coeffs = rem.get(j, {})
            # x_a^j * coeffs -> quotient term x_a^{j-1} * coeffs ; subtract (x_a - y) x_a^{j-1} coeffs
            for kk, cf in coeffs.items():
                if cf == 0: continue
                qk = list(kk); qk[a] = j - 1; q[tuple(qk)] = q.get(tuple(qk), 0) + cf
                # add y * x_a^{j-1} * cf to rem[j-1]
                yk = list(kk)
                if b != -1: yk[b] += 1
                yk = tuple(yk)
                rem.setdefault(j - 1, {})
                rem[j - 1][yk] = rem[j - 1].get(yk, 0) + cf
            rem[j] = {}
        if any(v != 0 for v in rem.get(low, {}).values()): return None
        return {k: v for k, v in q.items() if v != 0}
    while True:
        found = False
        for a in range(n):
            for b in [-1] + list(range(n)):
                if b == a: continue
                if vanishes_on(f, a, b):
                    q = divide(f, a, b)
                    if q is not None:
                        f = q; count += 1; found = True; break
            if found: break
        if not found: break
    return count, len(f) == 1

rng = random.Random(int(sys.argv[2]) if len(sys.argv) > 2 else 5)
trials = int(sys.argv[1]); ntight = 0; nonroot = []
for t in range(trials):
    n = rng.randint(1, int(sys.argv[3]) if len(sys.argv) > 3 else 3); R = rng.choice([1, 2, 3])
    S = sorted({tuple(rng.randint(-R, R) for _ in range(n)) for _ in range(rng.randint(2, int(sys.argv[4]) if len(sys.argv) > 4 else 14))})
    mw = MW(S, n)
    if mw != int(mw): continue
    m = int(mw)
    if m == 0: continue
    B = kernel_order(S, n, m)        # measures of order >= m = MW
    if not B: continue
    # generic element of the kernel
    coef = [Fr(rng.randint(-5, 5)) for _ in B]
    c = [sum(a * v[j] for a, v in zip(coef, B)) for j in range(len(S))]
    if all(x == 0 for x in c): continue
    supp = [S[j] for j in range(len(S)) if c[j] != 0]
    if MW(supp, n) != m: continue   # support shrank: then ord = m < MW(supp)? impossible by Thm B unless equal
    ntight += 1
    cnt, mono = strip_roots(c, S, n, rng)
    if not (mono and cnt == m):
        nonroot.append((S, c, cnt, mono))
print(f"tight measures found: {ntight}; not a root-binomial product: {len(nonroot)}")
for x in nonroot[:3]: print("  ", x)
