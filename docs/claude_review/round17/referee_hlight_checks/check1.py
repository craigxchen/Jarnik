import numpy as np, itertools, sys
from rlib import *
rng = np.random.default_rng(17)
out = []
def formula(q, ch, lam):
    f = lam[1:]; F = f.sum(); h = F + f
    c = np.array([sum(int(f[t])*int(ch[(a-t) % q]) for t in range(q)) for a in range(q)])
    return 2*(np.maximum(np.abs(h), np.abs(c)).sum() - q)
primes = [p for p in range(7, 90) if p % 4 == 3 and all(p % d for d in range(2, p))]
for q in primes:
    P = patterns(q); ch = chi_table(q); bad = 0
    # pair margins
    pm = set()
    for x in range(q+1):
        for y in range(x+1, q+1):
            lam = np.zeros(q+1, dtype=np.int64); lam[x] = 1; lam[y] = -1
            pm.add(sigma_def(P, lam))
    # AGL(1,q) invariance of pattern set (as partitions)
    canon = lambda r: tuple(r) if r[0] == 0 else tuple(1-r)
    Pset = {canon(r) for r in P}; inv = True
    for t in range(1, q):
        for b in range(q):
            perm = [0] + [1 + (t*x + b) % q for x in range(q)]
            for r in P:
                rr = np.zeros(q+1, dtype=np.int64); rr[perm] = r
                if canon(rr) not in Pset: inv = False
            if not inv: break
        if not inv: break
    for trial in range(60):
        lam = rng.integers(-3, 4, size=q+1); lam[0] -= lam.sum()
        if trial % 3 == 0:
            lam = np.zeros(q+1, dtype=np.int64); idx = rng.choice(q+1, 6, replace=False); lam[idx] = [1, 1, 1, -1, -1, -1]
        if sigma_def(P, lam) != formula(q, ch, lam): bad += 1
    out.append(f"q={q}: 2*sigma_pairs={sorted(pm)} AGL-invariant={inv} formula mismatches={bad}/60")
    print(out[-1], flush=True)
open('check1_out.txt', 'w').write('\n'.join(out) + '\n')
