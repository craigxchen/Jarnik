"""Adversarial search for near-flat characters in the medium range ||c||_1 in [M/4, 2.5M]
(the range that outside.md Sec. 2.2 left open).  Hill climbing on F(c) = ||H_A^T c||_1 restricted to
bands of n = ||c||_1, from structured starting points (restricted columns / half-sums, perturbed
columns / half-sums, random ternary, sums of columns).  Exceptional vectors (pairs, +-columns,
+-half-sums) are excluded.  Evidence only (Lemma S is proved); reports min (F - M)/M per band."""
import random
import numpy as np
from paley import paley, transform, exceptional

random.seed(77)
np.random.seed(77)


def starts(H, n_lo, n_hi, count):
    M = H.shape[0]
    out = []
    tries = 0
    while len(out) < count and tries < 50 * count:
        tries += 1
        t = random.randint(0, 5)
        c = np.zeros(M, dtype=np.int64)
        if t == 0:
            a = random.randint(1, M - 1)
            h = H[:, a]
            P = [x for x in range(M) if h[x] == 1]; N_ = [x for x in range(M) if h[x] == -1]
            k = random.randint(1, M // 2)
            c[random.sample(P, k)] = 1; c[random.sample(N_, k)] = -1
        elif t == 1:
            a, b = random.sample(range(1, M), 2)
            c = (H[:, a] + random.choice([1, -1]) * H[:, b]) // 2
            for _ in range(random.randint(1, 6)):
                x, y = random.sample(range(M), 2); c[x] += 1; c[y] -= 1
        elif t == 2:
            a = random.randint(1, M - 1)
            c = H[:, a].copy()
            for _ in range(random.randint(1, 6)):
                x, y = random.sample(range(M), 2); c[x] += 1; c[y] -= 1
        elif t == 3:
            lo_s, hi_s = max(2, n_lo), min(M, n_hi)
            if lo_s <= hi_s:
                s = random.randint(lo_s, hi_s) & ~1
                supp = random.sample(range(M), s)
                c[supp[: s // 2]] = 1; c[supp[s // 2:]] = -1
            else:  # large l1 norm: column plus random +-1 doubling on a subset
                a = random.randint(1, M - 1)
                c = H[:, a].copy()
                k = random.randint(1, M // 2)
                P = [x for x in range(M) if c[x] == 1]; N_ = [x for x in range(M) if c[x] == -1]
                c[random.sample(P, k)] += 1; c[random.sample(N_, k)] -= 1
        elif t == 4:
            k = random.randint(2, 3)
            cols = random.sample(range(1, M), k)
            v = sum(random.choice([1, -1]) * H[:, a] for a in cols)
            c = v // 2 if (v % 2 == 0).all() else v
        else:
            a, b = random.sample(range(1, M), 2)
            g = (H[:, a] + H[:, b]) // 2
            P = [x for x in range(M) if g[x] == 1]; N_ = [x for x in range(M) if g[x] == -1]
            k = random.randint(1, min(len(P), len(N_)))
            c[random.sample(P, k)] = 1; c[random.sample(N_, k)] = -1
            x, y = random.sample(range(M), 2); c[x] += 1; c[y] -= 1
        n = int(np.abs(c).sum())
        if c.any() and c.sum() == 0 and n_lo <= n <= n_hi and not exceptional(H, c):
            out.append(c)
    return out


def climb(H, c, n_lo, n_hi, steps=80):
    M = H.shape[0]
    Hn = H[:, 1:]
    D = (Hn[:, None, :] - Hn[None, :, :]).reshape(M * M, M - 1)
    T = transform(H, c)
    best = int(np.abs(T).sum())
    for _ in range(steps):
        Fs = np.abs(T[None, :] + D).sum(axis=1)
        order = np.argsort(Fs, kind="stable")
        moved = False
        for idx in order[:200]:
            if Fs[idx] >= best:
                break
            x, y = divmod(int(idx), M)
            if x == y:
                continue
            c2 = c.copy(); c2[x] += 1; c2[y] -= 1
            n2 = int(np.abs(c2).sum())
            if not (n_lo <= n2 <= n_hi) or not c2.any() or exceptional(H, c2):
                continue
            c, T, best = c2, T + D[idx], int(Fs[idx])
            moved = True
            break
        if not moved:
            break
    return c, best


def main():
    for q in (43, 47, 59, 67):
        H = paley(q)
        M = q + 1
        bands = [(M // 4, M // 3), (M // 3, int(0.45 * M)), (int(0.45 * M), int(0.55 * M)),
                 (int(0.55 * M), int(0.9 * M)), (int(0.9 * M), int(1.1 * M)),
                 (int(1.1 * M), int(1.6 * M)), (int(1.6 * M), int(2.5 * M))]
        line = []
        for lo, hi in bands:
            best = 10 ** 9
            for c in starts(H, lo, hi, 60):
                c2, f = climb(H, c, lo, hi)
                best = min(best, f - M)
            line.append("[%d,%d]: %.3f" % (lo, hi, best / M))
        print("q=%d M=%d  min (F-M)/M per n-band:  %s" % (q, M, "  ".join(line)))


if __name__ == "__main__":
    main()
