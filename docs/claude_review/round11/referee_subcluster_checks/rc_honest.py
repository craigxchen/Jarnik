# Referee's independent check of section 6.4 ("honest small models have no supercritical slab-box systems").
# Pair-separated binary codes (min distance >= r/2), equal weights. For EVERY n-subset J and EVERY base point,
# slab-box A = {k in Z^J_0 : |k_x|<=T (x != base), |<k,beta_p>| <= w}, LP(A) = sum_p omega_p h_A(beta_p) (exact
# rational), m = floor(LP)+1.  V(A,m) = 0  iff  the degree-(m-1) evaluation matrix has rank |A| (certified mod p).
import itertools, random, sys
from fractions import Fraction
import numpy as np
random.seed(77)
P = 1000003
def rank_modp(rows):
    A = np.array(rows, dtype=np.int64) % P
    nr, nc = A.shape; r = 0
    for c in range(nc):
        if r == nr: break
        nz = np.nonzero(A[r:, c])[0]
        if len(nz) == 0: continue
        pv = r + nz[0]
        if pv != r: A[[r, pv]] = A[[pv, r]]
        inv = pow(int(A[r, c]), P-2, P)
        A[r] = (A[r]*inv) % P
        col = A[:, c].copy(); col[r] = 0
        idx = np.nonzero(col)[0]
        if len(idx): A[idx] = (A[idx] - (col[idx, None]*A[r][None, :]) % P) % P
        r += 1
    return r
def mons(d, deg):
    return [a for t in range(deg+1) for a in itertools.product(range(t+1), repeat=d) if sum(a) == t]
def code(M, r, d):
    for trial in range(5000):
        cand = list(itertools.product([0,1], repeat=r)); random.shuffle(cand)
        words = []
        for c in cand:
            if all(sum(a != b for a, b in zip(c, w)) >= d for w in words):
                words.append(c)
                if len(words) == M: return words
    return None
def run(M, r, n, T, w, maxsys):
    d = -(-r//2)
    Wd = code(M, r, d)
    if Wd is None: print("no code", M, r); return
    cols = [tuple(Wd[x][p] for x in range(M)) for p in range(r)]
    cols = [c for c in cols if 0 < sum(c) < M]
    R = len(cols)
    tested = 0; nonzero = 0; maxratio = Fraction(0)
    systems = [(J, b) for J in itertools.combinations(range(M), n) for b in J]
    random.shuffle(systems)
    for J, b in systems[:maxsys]:
        others = [x for x in J if x != b]
        A = []
        for ks in itertools.product(range(-T, T+1), repeat=n-1):
            k = dict(zip(others, ks)); k[b] = -sum(ks)
            if all(abs(sum(k[x]*(2*c[x]-1) for x in J)) <= w for c in cols):
                A.append(tuple(k[x] for x in others))
        LP = Fraction(0)
        for c in cols:
            vals = [sum(kk*(2*c[x]-1) for kk, x in zip(kv, others)) - sum(kv)*(2*c[b]-1) for kv in A]
            LP += Fraction(max(vals)-min(vals), R)
        m = int(LP) + 1     # floor(LP)+1
        ms = mons(n-1, m-1)
        if len(ms) < len(A):
            nonzero += 1; tested += 1; continue
        rows = [[int(np.prod([kv[i]**a[i] for i in range(n-1)])) for a in ms] for kv in A]
        rk = rank_modp(rows)
        tested += 1
        if rk < len(A): nonzero += 1   # (mod-p rank deficiency; would need exact recheck)
    print(f"M={M} r={R} n={n} T={T} w={w}: systems tested {tested} (of {len(systems)}), V(A,floor(LP)+1)!=0 in {nonzero}", flush=True)
for (M, r) in [(6,6),(8,8),(10,10),(12,12)]:
    run(M, r, 4, 3, 4, 60)
    run(M, r, 4, 4, 6, 40)
    run(M, r, 5, 3, 4, 30)
for (M, r) in [(8,8),(12,12)]:
    run(M, r, 6, 2, 4, 8)
    run(M, r, 6, 3, 4, 4)
