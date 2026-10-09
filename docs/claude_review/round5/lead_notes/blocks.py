# Multi-block random-permutation profiles: per-block excess E(c) = ||H'^T Pi^T c||_1 - (M-1)
# for small-support characters c (sum c = 0). Exact integer arithmetic.
import numpy as np, itertools, random, sys
def sylvester(m):
    H = np.array([[1]])
    for _ in range(m): H = np.block([[H, H], [H, -H]])
    return H
def paley(q):  # q = 3 mod 4 prime, order q+1
    chi = [0]*q
    for x in range(1, q): chi[(x*x) % q] = 1
    chi = [1 if c else -1 for c in chi]; chi[0] = 0
    Q = np.array([[chi[(j-i) % q] if i != j else 0 for j in range(q)] for i in range(q)])
    S = np.zeros((q+1, q+1), dtype=int); S[0,1:] = 1; S[1:,0] = -1; S[1:,1:] = Q
    H = S + np.eye(q+1, dtype=int)
    # normalise first column to all ones
    H = H * H[:,[0]]
    return H
def excess(H, c):
    T = H.T @ c
    return int(np.abs(T[1:]).sum()) - (H.shape[0]-1)
def survey(H, name, trials=20000, supports=(2,3,4,5,6,8), seed=1):
    M = H.shape[0]; rng = random.Random(seed)
    print(f"== {name} M={M}")
    for s in supports:
        hist = {}
        for _ in range(trials):
            rows = rng.sample(range(M), s)
            # +-1 character with sum 0 (s even) or one entry 2 (s odd)
            c = np.zeros(M, dtype=int)
            if s % 2 == 0:
                for k, x in enumerate(rows): c[x] = 1 if k < s//2 else -1
            else:
                c[rows[0]] = 2
                for k, x in enumerate(rows[1:]): c[x] = 1 if k < (s-1)//2 - 1 + 0 else -1
                # fix sum: make sum zero
                tot = c.sum()
                if tot != 0:
                    # adjust: entries 2, then (s-1) entries: need sum(rest) = -2
                    c[:] = 0; c[rows[0]] = 2
                    neg = (s - 1 + 2)//2; pos = s - 1 - neg
                    for k, x in enumerate(rows[1:]): c[x] = -1 if k < neg else 1
            assert c.sum() == 0
            e = excess(H, c)
            hist[e] = hist.get(e, 0) + 1
        items = sorted(hist.items())
        print(f" s={s} n={int(np.abs(c).sum())}: min excess {items[0][0]} (freq {items[0][1]/trials:.4f}); "
              f"P(excess<=M/8)={sum(v for k,v in items if k <= M//8)/trials:.4f}; median {sorted(sum([[k]*v for k,v in items],[]))[trials//2]}")
if __name__ == "__main__":
    for m in (5, 6):
        survey(sylvester(m), f"Sylvester 2^{m}")
    for q in (31, 59):
        survey(paley(q), f"Paley q={q}")
