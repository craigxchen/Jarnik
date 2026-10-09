"""Theorem 1 on large Paley orders (exact integer arithmetic; evidence complementing the proof):
all pairs, all +-columns, all two-column half-sums (H_a +- H_b)/2, random +-1 characters of every support size,
and adversarial hill-climbing on  G(c) - bound(c)  over +-1 characters.  Also the same families for a RANDOM
row->label assignment (no skew pairing), to show what the pairing buys on rank-one rectangles."""
import sys, time
import numpy as np
from lib import *

def bound(c, b):
    n = int(np.abs(c).sum()); K = int(np.abs(c).max())
    if K >= 2:
        return (b - 3) * n / 2 + b
    return b - 4 if n == 2 else 2 * n + b - 8

def run(q, b, assign_kind, seed=0, climbs=40, steps=3000):
    H = normalize(paley1(q)); M = q + 1
    rng = np.random.default_rng(seed)
    if assign_kind == "skew":
        assign = np.array([1] + list(range(1, q + 1)))
    else:
        assign = np.concatenate([np.arange(1, M), [rng.integers(1, M)]]); rng.shuffle(assign)
    S, lab, frow = flipped_profile(H, b, assign)
    St = S.T.astype(np.int64); r = S.shape[1]
    G = lambda c: int(np.abs(St @ c).sum()) - r
    worst = (1e9, None, None)
    def upd(c, tag):
        nonlocal worst
        s = G(c) - bound(c, b)
        if s < worst[0]:
            worst = (s, tag, int(np.abs(c).sum()))
    # pairs
    for x in range(M):
        for y in range(x + 1, M):
            c = np.zeros(M, dtype=np.int64); c[x] = 1; c[y] = -1; upd(c, "pair")
    # columns and half-sums
    for a in range(1, M):
        upd(H[:, a].copy(), "column")
        for bb in range(a + 1, M):
            for s in (1, -1):
                c = (H[:, a] + s * H[:, bb]) // 2; upd(c, "half-sum")
    # random +-1
    for _ in range(20000):
        n = 2 * rng.integers(1, M // 2 + 1)
        idx = rng.choice(M, size=n, replace=False)
        c = np.zeros(M, dtype=np.int64); c[idx[: n // 2]] = 1; c[idx[n // 2:]] = -1; upd(c, "random")
    # adversarial hill-climb (moves: swap a +1/-1 position, move a +-1 to an empty position, add/remove a +-pair)
    for t in range(climbs):
        n = 2 * rng.integers(2, M // 2 + 1)
        idx = rng.choice(M, size=n, replace=False)
        c = np.zeros(M, dtype=np.int64); c[idx[: n // 2]] = 1; c[idx[n // 2:]] = -1
        if t % 4 == 0:   # start near a column or half-sum
            a = rng.integers(1, M); c = H[:, a].copy()
        cur = G(c) - bound(c, b)
        for _ in range(steps):
            d = c.copy()
            mv = rng.integers(0, 3)
            if mv == 0:
                x = rng.integers(0, M); y = rng.integers(0, M); d[x], d[y] = d[y], d[x]
            elif mv == 1:
                nz = np.flatnonzero(d); z = np.flatnonzero(d == 0)
                if len(z) == 0: continue
                x = rng.choice(nz); y = rng.choice(z); d[y] = d[x]; d[x] = 0
            else:
                z = np.flatnonzero(d == 0)
                if len(z) >= 2 and rng.random() < 0.5:
                    x, y = rng.choice(z, size=2, replace=False); d[x] = 1; d[y] = -1
                else:
                    p = np.flatnonzero(d == 1); m = np.flatnonzero(d == -1)
                    if len(p) and len(m) and len(p) + len(m) > 2:
                        d[rng.choice(p)] = 0; d[rng.choice(m)] = 0
            if not d.any(): continue
            val = G(d) - bound(d, b)
            if val <= cur:
                c, cur = d, val
        upd(c, "climb")
    return worst

if __name__ == "__main__":
    t0 = time.time()
    for q, b in [(43, 8), (103, 8), (199, 8), (43, 5), (199, 5)]:
        w1 = run(q, b, "skew")
        w2 = run(q, b, "random", seed=3)
        print("q=%3d M=%3d b=%d  skew pairing: min[G - bound] = %s (%s, n=%s)   random assignment: %s (%s, n=%s)   [%.0fs]"
              % (q, q + 1, b, w1[0], w1[1], w1[2], w2[0], w2[1], w2[2], time.time() - t0), flush=True)
