"""Search for instances with max block-zonotope order < m (no product certificate at 0) but dim V(X,m) > 0,
and record whether SL (char0 / p-integral) holds.  Blocks up to size 4 searched exhaustively."""
import itertools, random, sys
import numpy as np
from bool_sl import build, colmat, rank_modP, member_Zp, subsets_upto

def max_blocks(M, U, c, s, cap, bmax=4):
    cand = []
    for k in range(1, bmax + 1):
        cand += list(itertools.combinations(range(M), k))
    r = len(U)
    W = {B: np.array([sum(u[i] for i in B) for u in U]) for B in cand}
    cvec = np.array(c)
    # a single block must itself satisfy constraint
    cand = [B for B in cand if np.all(np.maximum(W[B], 0) <= s + cvec) and np.all(np.minimum(W[B], 0) >= cvec - s)]
    order = sorted(cand, key=lambda B: (np.abs(W[B]).sum(), len(B)))
    best = 0
    def dfs(start, used, pos, neg, L):
        nonlocal best
        if L > best: best = L
        if best >= cap: return
        if L + (M - sum(used)) <= best: return
        for idx in range(start, len(order)):
            B = order[idx]
            if any(used[i] for i in B): continue
            w = W[B]
            npos = pos + np.maximum(w, 0); nneg = neg + np.minimum(w, 0)
            if np.all(npos <= s + cvec) and np.all(nneg >= cvec - s):
                for i in B: used[i] = True
                dfs(idx + 1, used, npos, nneg, L + 1)
                for i in B: used[i] = False
                if best >= cap: return
    dfs(0, [False] * M, np.zeros(r), np.zeros(r), 0)
    return best

if __name__ == '__main__':
    rng = random.Random(int(sys.argv[1]))
    ntr = int(sys.argv[2])
    found = 0
    for trial in range(ntr):
        M = rng.choice([10, 11, 12])
        m = rng.choice([4, 5])
        r = rng.choice([4, 5, 6, 7, 8])
        s = 1
        U = [np.array([rng.choice([-1, 0, 1]) for _ in range(M)]) for _ in range(r)]
        c = [0] * r
        L = max_blocks(M, U, c, s, m)
        if L >= m:
            continue
        pts = build(M, m, U, c, s)
        Js = subsets_upto(M, m)
        zero = tuple([0] * M); pts.remove(zero); pts = [zero] + pts
        G = colmat(pts, Js)
        r_all = rank_modP(G)
        if r_all == len(pts):
            print(dict(M=M, m=m, r=r, X=len(pts), dimV=0, maxblocks=L)); sys.stdout.flush()
            continue
        r_oth = rank_modP(G[:, 1:])
        char0 = (r_all == r_oth)
        pint = None
        if char0:
            pint, mv, _ = member_Zp(G[:, 1:], G[:, 0], 3, 19)
        print(dict(M=M, m=m, r=r, X=len(pts), E=len(Js), dimV=len(pts) - r_all, char0=char0, pint=pint, maxblocks=L)); sys.stdout.flush()
