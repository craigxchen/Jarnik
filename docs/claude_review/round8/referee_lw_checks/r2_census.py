"""Referee check R2/R3 (independent method): exact minimum and census of s(c) for the skew double core.

R2 (rigorous minimum).  s(c) <= M/2 for some nonzero zero-sum c iff some set Z of 2q - M/2 types
     together with the all-ones vector has rank < M over Q.  We test every such Z: rank mod a prime
     p = 2^31 - 1 equal to M implies rank M over Q.  (q = 7: C(14,10) sets; q = 11: C(22,16) sets.)
R3 (census, independent of lw_checks).  For each (M-2)-subset of types, the kernel of the
     (M-1) x M matrix [types; ones] is spanned by the vector of signed maximal minors.  Minors are
     computed in floating point (|minor| <= 11^5.5 < 6e5, so rounding is safe) and every rounded vector
     is then verified EXACTLY (integer A c = 0); s(c) is computed exactly.  Degenerate subsets (all
     minors 0) are counted: if there are none, every character with zero set of size >= M-2 is found.
Usage: python3 r2_census.py
"""
import itertools, sys, time
import numpy as np
from math import gcd
from r1_dcore import build

P = 2**31 - 1
out = []
def log(s):
    print(s); out.append(s)

def rank_modp(rows, ncols):
    A = [[x % P for x in r] for r in rows]
    rk = 0
    for col in range(ncols):
        piv = next((i for i in range(rk, len(A)) if A[i][col]), None)
        if piv is None: continue
        A[rk], A[piv] = A[piv], A[rk]
        inv = pow(A[rk][col], P-2, P)
        A[rk] = [(x*inv) % P for x in A[rk]]
        for i in range(len(A)):
            if i != rk and A[i][col]:
                f = A[i][col]
                A[i] = [(x - f*y) % P for x, y in zip(A[i], A[rk])]
        rk += 1
        if rk == ncols: break
    return rk

def r2(q):
    M = q + 1
    K, s1, s2 = build(q)
    types = s1 + s2
    n = len(types); size = n - M//2
    bad = 0; tot = 0
    for Z in itertools.combinations(range(n), size):
        tot += 1
        if rank_modp([types[j] for j in Z] + [[1]*M], M) < M:
            bad += 1
    log(f"R2 q={q} M={M}: {tot} sets of {size} types (= #types - M/2); sets with rank < M mod p: {bad}"
        f" -> {'no c with s(c) <= M/2, so min s >= M/2+1 (rigorous)' if bad == 0 else 'NOT settled'}")

def r3(q):
    M = q + 1
    K, s1, s2 = build(q)
    types = s1 + s2
    n = len(types)
    T = np.array(types, dtype=np.int64)
    found = {}
    degenerate = 0; failed = 0; tot = 0
    it = itertools.combinations(range(n), M - 2)
    B = 4000
    t0 = time.time()
    while True:
        chunk = list(itertools.islice(it, B))
        if not chunk: break
        tot += len(chunk)
        idx = np.array(chunk)
        A = np.concatenate([T[idx], np.ones((len(chunk), 1, M), dtype=np.int64)], axis=1)  # (b, M-1, M)
        Af = A.astype(float)
        minors = np.empty((len(chunk), M))
        for k in range(M):
            sub = np.delete(Af, k, axis=2)
            minors[:, k] = ((-1)**k) * np.linalg.det(sub)
        C = np.rint(minors).astype(np.int64)
        resid = np.einsum('bij,bj->bi', A, C)
        zero = ~C.any(axis=1)
        degenerate += int(zero.sum())
        okv = (~resid.any(axis=1)) & (~zero)
        failed += int(((resid.any(axis=1)) & (~zero)).sum())
        img = C[okv] @ T.T
        svals = (img != 0).sum(axis=1)
        for cvec, s in zip(C[okv], svals):
            g = 0
            for x in cvec: g = gcd(g, abs(int(x)))
            cv = tuple(int(x)//g for x in cvec)
            if next(x for x in cv if x != 0) < 0: cv = tuple(-x for x in cv)
            found[cv] = int(s)
    # exact re-verification of every stored c (pure Python integers)
    for cv, s in found.items():
        assert sum(cv) == 0
        assert sum(1 for t in types if sum(a*b for a, b in zip(cv, t)) != 0) == s
    hist = {}
    for s in found.values(): hist[s] = hist.get(s, 0) + 1
    smin = min(hist)
    light = [cv for cv, s in found.items() if s == smin]
    cols = set()
    for fam in (s1, s2):
        for v in fam:
            cv = tuple(v)
            if next(x for x in cv if x != 0) < 0: cv = tuple(-x for x in cv)
            cols.add(cv)
    log(f"R3 q={q} M={M}: {tot} subsets of size M-2, degenerate (all minors 0): {degenerate}, "
        f"rounding failures: {failed}, distinct primitive characters: {len(found)} ({time.time()-t0:.0f}s)")
    log(f"   histogram s -> count (s <= M/2+5): {[(s, hist[s]) for s in sorted(hist) if s <= M//2 + 5]}")
    log(f"   min s = {smin}; lightest characters = the 2q type vectors themselves: {set(light) == cols}")

if __name__ == '__main__':
    for q in (7, 11):
        r2(q)
    for q in (7, 11):
        r3(q)
    open('r2_census_out.txt', 'w').write('\n'.join(out) + '\n')
