"""Boolean SL evidence, with product-certificate control.
For each instance also compute the largest 'block zonotope' at 0 inside X: disjoint nonempty blocks B_1..B_L
(v_i = 1_{B_i}) such that every union of blocks lies in X.  A product certificate of order m at 0 exists iff L >= m.
We report instances where char0/pint SL holds but max L < m (genuinely non-product certificates)."""
import itertools, random, sys
import numpy as np
from bool_sl import build, colmat, rank_modP, member_Zp, subsets_upto

def in_X(t, U, c, s):
    return all(abs(int(np.dot(t, u)) - cq) <= s for u, cq in zip(U, c))

def max_blocks(M, U, c, s, cap):
    # exhaustive search over block systems is expensive; we search singletons + pairs greedily-exhaustively:
    # blocks of size <= 2.  Union of blocks in X for all subsets: since the constraint is |sum over chosen blocks
    # of w_b| <= s with w_b = (<1_B,u_q>)_q, need: for each q, sum of positive w_b,q <= s+c_q and sum of
    # negative >= c_q - s.
    cand = [(i,) for i in range(M)] + list(itertools.combinations(range(M), 2))
    W = {B: np.array([sum(u[i] for i in B) for u in U]) for B in cand}
    best = 0
    # DFS over disjoint blocks maintaining positive/negative sums per q
    r = len(U)
    order = sorted(cand, key=lambda B: np.abs(W[B]).sum())
    def dfs(start, used, pos, neg, L):
        nonlocal best
        best = max(best, L)
        if best >= cap:
            return
        for idx in range(start, len(order)):
            B = order[idx]
            if any(used[i] for i in B):
                continue
            w = W[B]
            npos = pos + np.maximum(w, 0); nneg = neg + np.minimum(w, 0)
            if np.all(npos <= s + np.array(c)) and np.all(nneg >= np.array(c) - s):
                for i in B: used[i] = True
                dfs(idx + 1, used, npos, nneg, L + 1)
                for i in B: used[i] = False
                if best >= cap:
                    return
    dfs(0, [False] * M, np.zeros(r), np.zeros(r), 0)
    return best

if __name__ == '__main__':
    rng = random.Random(int(sys.argv[1]) if len(sys.argv) > 1 else 0)
    stats = []
    for trial in range(int(sys.argv[2]) if len(sys.argv) > 2 else 30):
        M = rng.choice([9, 10, 11, 12])
        m = rng.choice([4, 5])
        r = rng.choice([3, 4, 5, 6, 8])
        s = 1
        U = [np.array([rng.choice([-1, 0, 1]) for _ in range(M)]) for _ in range(r)]
        c = [0] * r
        pts = build(M, m, U, c, s)
        Js = subsets_upto(M, m)
        if len(pts) <= len(Js) // 2:
            continue
        zero = tuple([0] * M); pts.remove(zero); pts = [zero] + pts
        G = colmat(pts, Js)
        r_all = rank_modP(G); r_oth = rank_modP(G[:, 1:])
        char0 = (r_all == r_oth)
        L = max_blocks(M, U, c, s, m)
        pint = None
        if char0:
            pint, mv, _ = member_Zp(G[:, 1:], G[:, 0], 3, 19)
        line = dict(M=M, m=m, r=r, X=len(pts), E=len(Js), rank=r_all, dimV=len(pts) - r_all, char0=char0, pint=pint, maxblocks=L)
        print(line); sys.stdout.flush()
