"""Referee check: vertical polynomial templates of degree D as CIRCLE configurations (weakest.md Sec. 7.3).

Points z_sigma(t) = prod_l (t + i b_l)^{(1+sigma_l)/2} conj(t + i b_l)^{(1-sigma_l)/2}, sigma in {+-1}^D,
all of norm prod_l (t^2 + b_l^2).  arg z_sigma = sum_l sigma_l arctan(b_l/t)
= sum_{r odd} (-1)^((r-1)/2) (sum_l sigma_l b_l^r) / (r t^r).
A fibre {sigma : sum_l sigma_l b_l^r = c_r for all odd r < D/2} lies on an arc of length O(sqrt R)
(R = sqrt N ~ t^(D/2)), i.e. bounded C, for every large t.  With private primes (t in a suitable
residue class) the rank of the configuration is the affine rank of the fibre's sign vectors.

For each D and each b-subset of {1..B} we report the largest fibre and the largest affine rank of a
fibre (exact integer arithmetic; rank over Q via fractions).  Usage: python3 vertical_fibres.py D B
"""
import sys
import itertools
from fractions import Fraction
import numpy as np


def affine_rank(vecs):
    v0 = vecs[0]
    rows = [[Fraction(a - b) for a, b in zip(v, v0)] for v in vecs[1:]]
    if not rows:
        return 0
    m = rows
    rk, ncol = 0, len(m[0])
    for c in range(ncol):
        piv = next((i for i in range(rk, len(m)) if m[i][c] != 0), None)
        if piv is None:
            continue
        m[rk], m[piv] = m[piv], m[rk]
        for i in range(len(m)):
            if i != rk and m[i][c] != 0:
                f = m[i][c] / m[rk][c]
                m[i] = [a - f * b for a, b in zip(m[i], m[rk])]
        rk += 1
    return rk


def main(D, B, batch=4096):
    odd_r = [r for r in range(1, (D + 1) // 2 + 1) if r % 2 == 1 and 2 * r < D]
    S = np.array(list(itertools.product((-1, 1), repeat=D)), dtype=np.int64)  # 2^D x D
    best_size = (0, None, None)
    best_rank = (0, None, None, 0)
    hist = {}
    combos = itertools.combinations(range(1, B + 1), D)
    done = 0
    while True:
        chunk = list(itertools.islice(combos, batch))
        if not chunk:
            break
        bs = np.array(chunk, dtype=np.int64)  # nb x D
        key = np.zeros((S.shape[0], bs.shape[0]), dtype=np.int64)
        for r in odd_r:
            vals = S @ (bs ** r).T  # 2^D x nb
            span = 2 * int((bs ** r).sum(axis=1).max()) + 1
            key = key * span + (vals + span // 2)
        ks = np.sort(key, axis=0)
        # run lengths per column
        for j in range(bs.shape[0]):
            col = ks[:, j]
            change = np.flatnonzero(np.diff(col)) + 1
            starts = np.concatenate(([0], change))
            lens = np.diff(np.concatenate((starts, [len(col)])))
            mx = int(lens.max())
            hist[mx] = hist.get(mx, 0) + 1
            if mx >= 4:
                b = chunk[j]
                # examine every fibre of size >= 4 (rank)
                for st, ln in zip(starts, lens):
                    if ln < 4:
                        continue
                    kval = col[st]
                    idx = np.flatnonzero(key[:, j] == kval)
                    vecs = [list(S[i]) for i in idx]
                    rk = affine_rank(vecs)
                    if ln > best_size[0]:
                        best_size = (int(ln), b, rk)
                    if rk > best_rank[0] or (rk == best_rank[0] and ln > best_rank[3]):
                        best_rank = (rk, b, vecs, int(ln))
        done += len(chunk)
    print(f"D={D}, b subsets of {{1..{B}}} ({done} sets), conditions on odd r in {odd_r}")
    print("  histogram of max fibre size:", dict(sorted(hist.items())))
    print(f"  largest fibre: size {best_size[0]} at b={best_size[1]} (affine rank {best_size[2]})")
    print(f"  largest affine rank of a fibre (size>=4): {best_rank[0]} at b={best_rank[1]}, fibre size {best_rank[3]}")
    if best_rank[2]:
        for v in best_rank[2]:
            print("     ", v)


if __name__ == '__main__':
    main(int(sys.argv[1]), int(sys.argv[2]))
