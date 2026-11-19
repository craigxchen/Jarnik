"""Exact clusters for y^2 = 2t^2+1 over F3(t), with inert infinity.

For one pair of common high-coefficient strings, try all low parts of A,B
and sort their exact polynomial norms. Run: python3 vii_ramified_check.py 8
NumPy performs the batched lookup and sort with integer arithmetic.
"""
from argparse import ArgumentParser

import numpy as np


def coefficients(value, length):
    result = []
    for _ in range(length):
        value, digit = divmod(value, 3)
        result.append(digit)
    return result


def masks(poly):
    ones = sum(1 << j for j, c in enumerate(poly) if c % 3 == 1)
    twos = sum(1 << j for j, c in enumerate(poly) if c % 3 == 2)
    return ones, twos


def square_table(length, b_part):
    data = []
    for value in range(3**length):
        c = coefficients(value, length)
        square = [0] * (2*length-1)
        for j, cj in enumerate(c):
            for k, ck in enumerate(c):
                square[j+k] += cj*ck
        if b_part:
            # -d B^2=(t^2-1)B^2.
            square = [v-w for v, w in zip([0, 0]+square, square+[0, 0])]
        data.append(masks(square))
    return np.asarray(data, dtype=np.uint64)


def norm_keys(a, b, width):
    a1, a2 = a[:, 0, None], a[:, 1, None]
    b1, b2 = b[None, :, 0], b[None, :, 1]
    a0, b0 = ~(a1 | a2), ~(b1 | b2)
    ones = (a1 & b0) | (a0 & b1) | (a2 & b2)
    twos = (a2 & b0) | (a0 & b2) | (a1 & b1)
    return (ones | (twos << width)).ravel()


def trim(f):
    while f and f[-1] == 0:
        f.pop()
    return f


def is_squarefree(key, width):
    f = [((key >> j) & 1) + 2*((key >> (width+j)) & 1) for j in range(width)]
    g = trim([j*f[j] % 3 for j in range(1, width)])
    while g:
        while len(f) >= len(g):
            # Every nonzero residue modulo 3 is its own inverse.
            c = f[-1]*g[-1] % 3
            shift = len(f)-len(g)
            for j, gj in enumerate(g):
                f[shift+j] = (f[shift+j] - c*gj) % 3
            trim(f)
        f, g = g, f
    return len(f) == 1


def maxima(D):
    k = D//2
    width = 2*D+1
    if 2*width > 64:
        raise ValueError("The integer key has room for 2D+1 <= 32 coefficients.")
    sa = square_table(D+1, False)
    sb = square_table(D, True)
    low_a = 3 ** (k+1)
    low_b = 3 ** k
    overall = [0, 0, 0]
    squarefree = [0, 0, 0]
    for side, (a_lead, b_lead) in enumerate(((1, 0), (0, 1), (1, 1))):
        for a_start in range(a_lead*3**D, (a_lead+1)*3**D, low_a):
            a = sa[a_start:a_start+low_a]
            for b_start in range(b_lead*3**(D-1), (b_lead+1)*3**(D-1), low_b):
                b = sb[b_start:b_start+low_b]
                keys = norm_keys(a, b, width)
                values, counts = np.unique(keys, return_counts=True)
                best = int(counts.max())
                overall[side] = max(overall[side], best)
                if best > squarefree[side]:
                    for j in np.flatnonzero(counts > squarefree[side]):
                        if counts[j] > squarefree[side] and is_squarefree(int(values[j]), width):
                            squarefree[side] = int(counts[j])
    return overall, squarefree


if __name__ == "__main__":
    parser = ArgumentParser(description=__doc__)
    parser.add_argument("max_degree", nargs="?", default=8, type=int)
    args = parser.parse_args()
    print("D all monic parent_lc_2 squarefree squarefree_monic | (A,0),(0,B),(A,B)")
    for D in range(1, args.max_degree+1):
        overall, squarefree = maxima(D)
        print(D, max(overall), max(overall[:2]), overall[2], max(squarefree),
              max(squarefree[:2]), "|", overall, squarefree, flush=True)
