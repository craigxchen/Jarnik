#!/usr/bin/env python3
"""Supplement to ref_nontwin.py (1): a 0/1 profile (every column nonconstant, distinct rows) where the
literal recipe kappa' = d_M u (u injective mod n/g) gives a unit-1 pair collision of kappa = A nut + 2kappa'."""
import itertools, math, random
import importlib.util, sys
spec = importlib.util.spec_from_file_location('rn', 'ref_nontwin.py')
src = open('ref_nontwin.py').read().split('# ---------------- (1)')[0]
ns = {}; exec(src, ns); snf_diag = ns['snf_diag']
random.seed(9)
for trial in range(200000):
    M, r = random.choice([3, 4]), random.randint(4, 7)
    A = [[random.randint(0, 1) for _ in range(r)] for _ in range(M)]
    if len(set(map(tuple, A))) < M or any(len(set(A[i][j] for i in range(M))) == 1 for j in range(r)): continue
    d = snf_diag(A)
    if len(d) < M or d[-1] == 1: continue
    dM = d[-1]
    for n in (4, 6, 8, 10, 12):
        g = math.gcd(dM, n)
        if g == 1 or n // g < M: continue
        for nut in [tuple(random.randint(0,1) for _ in range(r)) for _ in range(12)]:
            base = [sum(A[x][j] * nut[j] for j in range(r)) for x in range(M)]
            for u in itertools.islice(itertools.permutations(range(n // g), M), 400):
                kp = [dM * ux % n for ux in u]
                kap = [(base[x] + 2 * kp[x]) % (2 * n) for x in range(M)]
                if len(set(kap)) < M:
                    print('0/1 profile A=%s (invariant factors %s), n=%d, nut=%s: kappa\'=%s injective mod n, kappa=%s mod %d collides'
                          % (A, d, n, nut, kp, kap, 2 * n))
                    sys.exit(0)
print('no 0/1 example found')
