"""Referee check: aligned Walsh flats on the Sylvester fake profiles (b = 6M), exact integers.

Aligned flats are the characters that break flipped-Walsh profiles (walsh-extraction Theorem C,
items 333/335): c_x = (-1)^{l(x - x0)} on a coset x0 + V of a subspace V of F_2^t, l a nonzero
functional on V.  They are flat: F(c) = M.  Every subspace V of F_2^t (t = 3, 4, 5), every coset and
every nonzero l is tested, with 5 independent random flip assignments per order, in the equal-weight
form of Lemma A.7(iii):  sum_j |c^T S_j| - r >= 2(n + M),  and F(c) >= M.
"""
import random
from itertools import product
import numpy as np
from r1_lemmaA7 import sylvester, normalize, build

random.seed(7)


def subspaces(t):
    """all subspaces of F_2^t, as sorted tuples of integer vectors (bitmasks)."""
    seen = set()
    out = []
    vecs = list(range(1, 2 ** t))
    # generate by closure of random-free enumeration over spanning sets of size <= t
    def span(gens):
        S = {0}
        for g in gens:
            S |= {s ^ g for s in S}
        return tuple(sorted(S))
    frontier = [()]
    for _ in range(t + 1):
        new = []
        for gens in frontier:
            sp = span(gens)
            if sp not in seen:
                seen.add(sp)
                out.append(sp)
                for v in vecs:
                    if v not in sp:
                        new.append(gens + (v,))
        # dedupe frontier by span
        d = {}
        for g in new:
            d.setdefault(span(g), g)
        frontier = list(d.values())
    return out


def main():
    for t in (3, 4, 5):
        M = 2 ** t
        H = normalize(sylvester(M))
        subs = [V for V in subspaces(t) if len(V) >= 2]
        worst = None
        count = 0
        for trial in range(5):
            S, cols, idx = build(H, 6 * M)
            r = S.shape[1]
            for V in subs:
                dim = len(V).bit_length() - 1
                # basis of V
                basis = []
                span = {0}
                for v in V:
                    if v not in span:
                        basis.append(v)
                        span |= {s ^ v for s in span}
                cosets = {}
                for x0 in range(M):
                    key = min(x0 ^ v for v in V)
                    cosets.setdefault(key, x0)
                for x0 in cosets.values():
                    for lvals in product((0, 1), repeat=dim):
                        if not any(lvals):
                            continue
                        c = np.zeros(M, dtype=np.int64)
                        for coeffs in product((0, 1), repeat=dim):
                            v = 0
                            for cb, b in zip(coeffs, basis):
                                if cb:
                                    v ^= b
                            par = sum(cb * lv for cb, lv in zip(coeffs, lvals)) % 2
                            c[x0 ^ v] = (-1) ** par
                        n = int(np.abs(c).sum())
                        T = c @ H
                        assert T[0] == 0
                        F = int(np.abs(T[1:]).sum())
                        assert F == M, F  # flat
                        G = int(np.abs(c @ S).sum()) - r
                        sl = G / 2 - (n + M)
                        assert sl >= 0
                        worst = sl if worst is None else min(worst, sl)
                        count += 1
        print("Sylvester M=%2d b=%3d: %d aligned flats (all subspaces, cosets, nonzero functionals, 5 random"
              " flip assignments): all flat (F=M) and Lemma A.7(iii) equal-weight form holds; min slack %.1f"
              % (M, 6 * M, count, worst))


if __name__ == "__main__":
    main()
