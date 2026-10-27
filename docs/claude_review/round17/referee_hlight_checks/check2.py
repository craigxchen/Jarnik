import numpy as np, itertools
from rlib import *
out = []
def log(s): print(s, flush=True); out.append(s)
def is_pair_mult(lam):
    nz = np.nonzero(lam)[0]
    return len(nz) == 2
# exhaustive: q=7 entries in [-2,2]; q=11 ternary.  All sum-zero lambda over q+1 coords (incl inf).
for q, B in [(7, 2), (7, 3), (11, 1)]:
    P = patterns(q); n = q+1
    vals = np.arange(-B, B+1)
    grid = np.array(list(itertools.product(vals, repeat=n-1)), dtype=np.int64)
    last = -grid.sum(1)
    keep = np.abs(last) <= B
    L = np.concatenate([grid[keep], last[keep][:, None]], axis=1)
    nzc = (L != 0).sum(1)
    S2 = 2*np.abs(L @ P.T).sum(1) - 2*q   # 2*sigma
    nonpair = nzc >= 3
    m = S2[nonpair].min()
    arg = L[nonpair][S2[nonpair] == m][0]
    log(f"q={q} |lam_x|<={B}: #lambda={len(L)}, min sigma over non-pair-multiples = {m/2} (r/6={2*q/6:.3f}), example {arg.tolist()}")
    # min sigma over multiples of pairs with |l|>=2
    mp = S2[(nzc == 2) & (np.abs(L).max(1) >= 2)]
    if len(mp): log(f"   multiples l*pair, |l|>=2: min sigma = {mp.min()/2} (claim |l|(q+1)-q -> {2*(q+1)-q})")
# base 4-terms e_x + e_y - e_z - e_inf: distribution of sigma
for q in [7, 11, 19, 23, 31, 43, 47]:
    P = patterns(q); n = q+1; vals = set()
    for x, y in itertools.combinations(range(1, n), 2):
        for z in range(1, n):
            if z in (x, y): continue
            lam = np.zeros(n, dtype=np.int64); lam[[x, y]] = 1; lam[z] = -1; lam[0] = -1
            vals.add(sigma_def(P, lam)/2)
    log(f"q={q}: base 4-terms e_x+e_y-e_z-e_inf sigma values = {sorted(vals)}  ((q-1)/2={(q-1)/2})")
open('check2_out.txt', 'w').write('\n'.join(out) + '\n')
