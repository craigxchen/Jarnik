import numpy as np, itertools
from rlib import *
out = []
def log(s): print(s, flush=True); out.append(s)
for q in [19, 23]:
    P = patterns(q); n = q+1
    # F=0, n=6, ternary on F_q ; |F|=1, n=6 (lam_inf=-1, f: three +1, two -1); n=8 F=0 sampled
    rows = []
    for S in itertools.combinations(range(1, n), 6):
        for plus in itertools.combinations(S[1:], 2):  # WLOG S[0] is +1 (negation symmetry)
            lam = np.zeros(n, dtype=np.int64); lam[list(S)] = -1; lam[S[0]] = 1; lam[list(plus)] = 1; rows.append(lam)
    L = np.array(rows); s = (np.abs(L @ P.T).sum(1) - q)
    log(f"q={q}: F=0, n=6 exhaustive ({len(L)}): min sigma={s.min()}  Prop2.5 bound (q-3)/2={(q-3)/2}, g(6)={(q*4+24-36)/8}")
    rows = []
    for S in itertools.combinations(range(1, n), 5):
        for minus in itertools.combinations(S, 2):
            lam = np.zeros(n, dtype=np.int64); lam[list(S)] = 1; lam[list(minus)] = -1; lam[0] = -1; rows.append(lam)
    L = np.array(rows); s = (np.abs(L @ P.T).sum(1) - q)
    log(f"q={q}: |F|=1, n=6 exhaustive ({len(L)}): min sigma={s.min()}  Prop2.4 bound (q+1)/2={(q+1)/2}")
    rows = []
    for x, y in itertools.combinations(range(1, n), 2):
        for z in range(1, n):
            if z in (x, y): continue
            lam = np.zeros(n, dtype=np.int64); lam[[x, y]] = 1; lam[z] = -1; lam[0] = -1; rows.append(lam)
    L = np.array(rows); s = (np.abs(L @ P.T).sum(1) - q)
    log(f"q={q}: |F|=1, n=4 exhaustive: min sigma={s.min()}  Prop2.4 bound (q-1)/2={(q-1)/2}")
# Prop 2.5(2) g(n) vs random F=0 ternary at q=43, n = 4..22
q = 43; P = patterns(q); n1 = q+1; rng = np.random.default_rng(5)
for nn in range(4, 23, 2):
    best = 10**9
    for _ in range(20000):
        lam = np.zeros(n1, dtype=np.int64); idx = rng.choice(np.arange(1, n1), nn, replace=False)
        lam[idx[:nn//2]] = 1; lam[idx[nn//2:]] = -1; best = min(best, int(np.abs(P @ lam).sum() - q))
    g = (q*(nn-2) + 4*nn - nn*nn)/(nn+2)
    log(f"q=43 F=0 n={nn}: sampled min sigma={best}; proved max(n-1, g(n) [n<=(q+1)/2]) = {max(nn-1, g if nn <= (q+1)//2 else -1):.2f}")
open('check5_out.txt', 'w').write('\n'.join(out) + '\n')
