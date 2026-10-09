"""Exact checks of Theorem 1 (character margin for identity-assigned fully flipped Paley profiles).

Profile P(q,b): H = normalised Paley-I matrix of order M = q+1 (row/column 0 = infinity, row/col i = element i-1
of F_q); b unflipped copies of every label a = 1..q; one flipped copy per row x, of label a(x) = x for
x = 1..q and a(0) = a_star.  G(c) = sum_j (|(S^T c)_j| - 1).
Claims checked (all integer arithmetic):
  (L) every 2x2 principal minor of H on rows/cols {i,j} subset {1..q} is nonzero (skew structure);
  (T1) G(c) >= b - 4 for every pair, G(c) >= 2n + b - 8 for every +-1 valued c with n = ||c||_1 >= 4;
  (T2) G(c) >= (b-3) n/2 + b for every c with ||c||_inf >= 2;
  (E)  the exact decomposition G = (b+1)(F-M) + b + |T_{a*}| + sum_x Delta_x.
Exhaustive over the stated ranges; the theorem itself is proved in construct.md.
"""
import sys, itertools, time
import numpy as np
from lib import *

def profile(q, b, a_star=1):
    H = normalize(paley1(q)); M = q + 1
    assign = np.array([a_star] + list(range(1, q + 1)))
    S, lab, frow = flipped_profile(H, b, assign)
    return H, S, assign

def check_minors(q):
    H = normalize(paley1(q))
    sub = H[1:, 1:]
    d = sub * sub.T  # H_ij H_ji
    off = d[~np.eye(q, dtype=bool)]
    diag = np.diag(sub)
    # det = H_ii H_jj - H_ij H_ji ; with diag entries all equal (=-1) and H_ij H_ji = -1 off-diagonal
    assert (off == -1).all() and (np.abs(diag) == 1).all()
    dets = np.outer(diag, diag) - d
    return int(np.abs(dets[~np.eye(q, dtype=bool)]).min())

def pm1_vectors(M, n):
    """all c in {-1,0,1}^M with sum 0 and exactly n nonzeros (n even), as an int8 array, in chunks"""
    half = n // 2
    signs = []
    for pos in itertools.combinations(range(n), half):
        s = -np.ones(n, dtype=np.int8); s[list(pos)] = 1; signs.append(s)
    signs = np.array(signs)
    buf = []
    for supp in itertools.combinations(range(M), n):
        c = np.zeros((len(signs), M), dtype=np.int8); c[:, supp] = signs
        buf.append(c)
        if len(buf) * len(signs) >= 200000:
            yield np.concatenate(buf); buf = []
    if buf:
        yield np.concatenate(buf)

def run(q, b, nmax):
    H, S, assign = profile(q, b)
    M = q + 1; r = S.shape[1]
    S16 = S.astype(np.int32)
    out = []
    for n in range(2, nmax + 1, 2):
        worst = None; cnt = 0
        for C in pm1_vectors(M, n):
            Y = C.astype(np.int32) @ S16
            G = np.abs(Y).sum(axis=1) - r
            bound = (b - 4) if n == 2 else (2 * n + b - 8)
            slack = G - bound
            i = int(np.argmin(slack)); cnt += len(C)
            if worst is None or slack[i] < worst[0]:
                worst = (int(slack[i]), int(G[i]), C[i].copy())
        assert worst[0] >= 0, (q, b, n, worst)
        Fw = F_of(H, worst[2].astype(np.int64))
        out.append((n, cnt, worst[1], worst[0], Fw))
    return out

def run_general(q, b, maxnorm, K=3):
    """all c with entries in [-K,K], sum 0, ||c||_1 <= maxnorm and ||c||_inf >= 2: check (T2)"""
    H, S, assign = profile(q, b)
    M = q + 1; r = S.shape[1]
    worst = None; cnt = 0
    vals = [v for v in range(-K, K + 1) if v != 0]
    for k in range(2, maxnorm + 1):
        for supp in itertools.combinations(range(M), k):
            for vv in itertools.product(vals, repeat=k):
                if sum(vv) != 0 or sum(abs(v) for v in vv) > maxnorm or max(abs(v) for v in vv) < 2:
                    continue
                c = np.zeros(M, dtype=np.int64); c[list(supp)] = vv
                n = int(np.abs(c).sum()); G = Gval(S, c); cnt += 1
                sl = G - ((b - 3) * n / 2 + b)
                assert sl >= 0, (q, b, c, G)
                if worst is None or sl < worst[0]:
                    worst = (sl, G, n)
    return cnt, worst

def check_decomposition(q, b, trials=2000, seed=1):
    H, S, assign = profile(q, b)
    M = q + 1; rng = np.random.default_rng(seed)
    a_star = assign[0]
    for _ in range(trials):
        k = rng.integers(2, M + 1)
        c = np.zeros(M, dtype=np.int64)
        idx = rng.choice(M, size=k, replace=False)
        c[idx] = rng.integers(-3, 4, size=k)
        c[idx[0]] -= c.sum()
        if not c.any():
            continue
        T = H.T @ c; F = int(np.abs(T[1:]).sum())
        h = np.array([H[x, assign[x]] for x in range(M)])
        Delta = np.abs(T[assign] - 2 * c * h) - np.abs(T[assign])
        rhs = (b + 1) * (F - M) + b + abs(int(T[a_star])) + int(Delta.sum())
        assert Gval(S, c) == rhs
    return trials

if __name__ == "__main__":
    t0 = time.time()
    print("(L) min |2x2 principal minor| over F_q x F_q:")
    for q in [7, 11, 19, 23, 31, 43, 47, 59, 67, 71, 79, 83, 103, 107, 127, 131, 139, 151, 163, 167, 179, 191, 199]:
        print("   q=%3d: %d" % (q, check_minors(q)))
    print("(E) exact decomposition on random integer c:")
    for q, b in [(11, 5), (19, 5), (23, 7), (43, 5), (67, 12)]:
        print("   q=%d b=%d: %d random c OK" % (q, b, check_decomposition(q, b)))
    print("(T1) exhaustive +-1 characters: n, #c, G at worst, slack over bound, F at worst")
    for q, b, nmax in [(11, 5, 12), (11, 8, 12), (19, 5, 8), (23, 5, 6), (31, 5, 6), (43, 5, 4)]:
        for row in run(q, b, nmax):
            print("   q=%d b=%d n=%2d  #c=%9d  G=%4d  slack=%3d  F=%d" % ((q, b) + row), flush=True)
    print("(T2) exhaustive c with entries in [-3,3], ||c||_inf >= 2, small norm:")
    for q, b, mx in [(11, 5, 8), (19, 5, 6), (23, 5, 6)]:
        cnt, w = run_general(q, b, mx)
        print("   q=%d b=%d ||c||_1<=%d: %d vectors, min slack %.1f (G=%d, n=%d)" % (q, b, mx, cnt, w[0], w[1], w[2]))
    print("total time %.0fs" % (time.time() - t0))
