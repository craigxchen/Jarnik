"""Check Theorem 4.1 of lw.md: for the Paley skew double core of order M = q+1, every nonzero
zero-sum integer c has s(c) >= M/2, where s(c) = number of the 2(M-1) cut types with <c,s_j> != 0.

Part A: the identities u_a + v_a = 2(c_a - c_0), u_a - v_a = 2 sum_y chi(a-y) c_y (exact, random c).
Part B: exhaustive minimum of s(c) over all zero-sum c with entries in {-1,0,1} (q = 7, 11) and
        with entries in [-2,2] and support <= 6 (q = 7, 11, 19).
Part C: exact maximum hyperplane load (= 2(M-1) - min s) by enumerating all hyperplanes spanned by
        M-2 types plus the all-ones vector (q = 7, 11); maximisers re-verified in exact arithmetic.
Usage: python3 check_dcore_support.py
"""
import itertools, random, sys
import numpy as np
from dcore import paley_skew, double_core_types, image, support_size, legendre, nullvec_int

random.seed(1)
out = []

def log(s):
    print(s)
    out.append(s)

# Part A
for q in (7, 11, 19, 23):
    K = paley_skew(q)
    M = q + 1
    f1, f2 = double_core_types(K)
    for trial in range(200):
        c = [random.randint(-5, 5) for _ in range(M - 1)]
        c.append(-sum(c))
        random.shuffle(c)
        u = image(c, f1)   # u_a for a = 1..q (cut vectors are K(x,a)K(x,0) = K(x,a))
        v = image(c, f2)
        for a in range(1, M):
            A = a - 1
            sig = sum(legendre(A - (y - 1), q) * c[y] for y in range(1, M))
            # cut of K^T has sign pattern K^T(x,a)K^T(x,0); its image differs from (K c~)_a by sign
            assert u[a - 1] == c[a] - c[0] + sig, (q, 'u')
            assert abs(v[a - 1]) == abs(c[a] - c[0] - sig), (q, 'v')
    log(f"Part A q={q}: identities u_a = (c_a-c_0)+Sigma_a, |v_a| = |(c_a-c_0)-Sigma_a| hold on 200 random c")

# Part B
def min_support_enum(q, entries, maxsupp):
    K = paley_skew(q)
    M = q + 1
    f1, f2 = double_core_types(K)
    T = np.array(f1 + f2, dtype=np.int64)  # (2q) x M
    best = None
    count = 0
    vals = [e for e in entries if e != 0]
    for s in range(2, maxsupp + 1):
        for supp in itertools.combinations(range(M), s):
            sub = T[:, supp]
            for signs in itertools.product(vals, repeat=s):
                if sum(signs) != 0:
                    continue
                if signs[0] < 0:
                    continue  # c and -c give the same support
                count += 1
                img = sub @ np.array(signs, dtype=np.int64)
                sc = int(np.count_nonzero(img))
                if best is None or sc < best[0]:
                    best = (sc, supp, signs)
    return best, count

for q, entries, maxsupp in ((7, (-1, 0, 1), 8), (11, (-1, 0, 1), 12), (7, (-2, -1, 0, 1, 2), 6),
                             (11, (-2, -1, 0, 1, 2), 6), (19, (-1, 0, 1), 6), (19, (-2, -1, 0, 1, 2), 4)):
    best, count = min_support_enum(q, entries, maxsupp)
    M = q + 1
    log(f"Part B q={q} M={M} entries={entries} supp<={maxsupp}: {count} characters (up to sign); "
        f"min s(c) = {best[0]} (bound M/2 = {M//2}); a minimiser: rows {best[1]} values {best[2]}")
    assert best[0] >= M // 2

# Part C
def max_load_exact(q):
    K = paley_skew(q)
    M = q + 1
    f1, f2 = double_core_types(K)
    types = f1 + f2
    T = np.array(types, dtype=float)
    n = len(types)
    ones = np.ones(M)
    best = -1
    maxis = []
    combs = list(itertools.combinations(range(n), M - 2))
    B = 20000
    for i0 in range(0, len(combs), B):
        chunk = np.array(combs[i0:i0 + B])
        A = np.concatenate([T[chunk], np.broadcast_to(ones, (len(chunk), 1, M))], axis=1)  # (b, M-1, M)
        U, S, Vt = np.linalg.svd(A)
        ok = S[:, -1] > 1e-8  # rank M-1
        c = Vt[:, -1, :]   # null vector
        img = np.abs(c @ T.T)
        loads = (img < 1e-7).sum(axis=1)
        loads[~ok] = -1
        m = loads.max()
        if m > best:
            best = m
            maxis = [tuple(chunk[j]) for j in np.nonzero(loads == m)[0][:50]]
        elif m == best:
            maxis += [tuple(chunk[j]) for j in np.nonzero(loads == m)[0][:50 - len(maxis)]] if len(maxis) < 50 else []
    # exact re-verification of maximisers
    exact = set()
    for sub in maxis:
        rows = [types[j] for j in sub] + [[1] * M]
        c = nullvec_int(rows, M)
        if c is None:
            continue
        load = sum(1 for t in image(c, types) if t == 0)
        exact.add((load, tuple(c)))
    return best, len(combs), sorted(exact, reverse=True)[:5]

for q in (7, 11):
    best, ncomb, ex = max_load_exact(q)
    M = q + 1
    log(f"Part C q={q} M={M}: {ncomb} hyperplanes enumerated; max load (float) = {best}; "
        f"min s(c) = {2*(M-1) - best} (bound M/2 = {M//2})")
    for load, c in ex[:3]:
        log(f"   exact maximiser: load {load}, s(c) = {2*(M-1)-load}, c = {list(c)}")

open('check_dcore_support_out.txt', 'w').write('\n'.join(out) + '\n')
