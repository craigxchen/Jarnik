"""Enumerate admissible model instances (n, profile, T, L, m):
 (i) slabs cut: A != B_T;  (ii) h_A(beta_p) < m for every column p;  (iii) V(A,m) != 0.
Uniform slab half-width L (the theorem uses one w for all columns: w = 2L).
Profiles = sets of sign patterns S (up to complement); full rank = the functionals k(S) span (R^d)^*.
"""
import itertools, sys, random
import numpy as np
from core import *

def patterns(n):
    pats = []
    for s in range(1, n // 2 + 1):
        for S in itertools.combinations(range(n), s):
            if 2 * s == n and 0 not in S:
                continue
            pats.append(S)
    return pats

def fullrank(n, prof):
    d = n - 1
    rows = []
    for S in prof:
        rows.append([ (1 if x in S else 0) - (1 if (n - 1) in S else 0) for x in range(d)])
    return np.linalg.matrix_rank(np.array(rows, dtype=float)) == d if rows else False

def admissible(n, prof, T, L, mmax=40):
    slabs = [(S, L) for S in prof]
    A = build_A(n, T, slabs)
    B = (2 * T + 1) ** (n - 1)
    if len(A) == B:
        return None
    hmax = max(2 * width(A, S) for S in prof)
    # Hilbert function: V(A,m) != 0 iff H_A(m-1) < |A|
    res = []
    m = hmax + 1
    while m <= mmax:
        r = moment_matrix(A, m).rank()
        if r == len(A):
            break
        res.append((m, len(A) - r))
        m += 1
    return dict(A=len(A), box=B, hmax=hmax, ms=res)

if __name__ == "__main__":
    n = int(sys.argv[1])
    pats = patterns(n)
    print("patterns", pats)
    out = []
    if n == 3 or n == 4:
        profs = [p for r in range(1, len(pats) + 1) for p in itertools.combinations(pats, r)]
    else:
        random.seed(1)
        profs = set()
        while len(profs) < int(sys.argv[2]) if len(sys.argv) > 2 else 40:
            r = random.randint(n, 3 * n)
            cols = [tuple(random.choice([1, -1]) for _ in range(n)) for _ in range(r)]
            prof = set()
            for c in cols:
                S = tuple(x for x in range(n) if c[x] == 1)
                Sc = tuple(x for x in range(n) if c[x] == -1)
                if len(S) == 0 or len(Sc) == 0:
                    continue
                key = min([S, Sc], key=lambda s: (len(s), [x != 0 for x in s] if 2*len(s)==n else 0, s))
                if 2 * len(key) == n and 0 not in key:
                    key = Sc if 0 in Sc else S
                prof.add(key)
            profs.add(tuple(sorted(prof)))
        profs = sorted(profs)
    Tmax = {3: 6, 4: 5, 5: 3}[n]
    for prof in profs:
        fr = fullrank(n, prof)
        for L in range(1, 4):
            for T in range(L, Tmax + 1):
                if (2 * T + 1) ** (n - 1) > 3000:
                    continue
                r = admissible(n, prof, T, L)
                if r and r['ms']:
                    print(prof, "fullrank" if fr else "deficient", "T=%d L=%d" % (T, L), r, flush=True)
