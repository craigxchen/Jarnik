# EXACT dims. Pair-separated squarefree full-rank profiles (binary codes, equal weights): do LP-supercritical
# sub-cluster systems (m > sum_p omega_p h_A(beta_p), V(A,m) != 0) exist at small n?  Also the pushforward to the
# global position space is the identity here (full rank), so 'combined' = sum of embedded V(A_J,m_J).
import itertools, random, math, sys
from engine import choose_m
random.seed(3)
def code(M, r, d):
    best = None
    for trial in range(20000):
        words = []
        cand = list(itertools.product([0, 1], repeat=r)); random.shuffle(cand)
        for c in cand:
            if all(sum(a != b for a, b in zip(c, w)) >= d for w in words): words.append(c)
            if len(words) == M: break
        if len(words) == M:
            return words
    return None
for (M, r) in [(6, 6), (8, 7), (8, 8), (10, 10), (12, 12)]:
    d = math.ceil(r/2)
    W = code(M, r, d)
    if W is None: print(f"M={M} r={r}: no code with distance >= {d}"); continue
    cols = [tuple(W[x][p] for x in range(M)) for p in range(r)]
    cols = [c for c in cols if 0 < sum(c) < M]
    print(f"M={M} r={len(cols)} pair-separated code (min dist {min(sum(a!=b for a,b in zip(u,v)) for u,v in itertools.combinations(W,2))} >= r/2)", flush=True)
    for n in [4, 5]:
        for (T, w) in [(3, 4), (4, 6), (4, 8)]:
            if n == 5 and T > 3: continue
            stats = []
            for J in list(itertools.combinations(range(M), n))[:6]:
                m, U, dV, lp, mh = choose_m(list(J), J[-1], cols, T, w, [1/len(cols)]*len(cols), 1)
                stats.append((U, lp, m, dV))
            print(f"   n={n} T={T} w={w}: (|A|, LP, m=floor(LP)+1, dim V(A,m)) for 6 subsets: {[(u, round(l,2), mm, dv) for u,l,mm,dv in stats]}", flush=True)
