#!/usr/bin/env python3
"""Own checks for the same-norm P2 comparison.
(a) normalisation: column sign counts on the class A, and survival of both exponents 0 and 2
    after deleting rows (worst case deterministic bound);
(b) Legendre parity = singleton indicator of the owner; Sidon property of unary signature sums;
(c) Reed-Muller dimension count;
(d) pair-margin statistic under a uniform permutation of values in an interval of width sigma:
    weighted sum sum_j c_j x_{pi(j)}, c_j = |K_ij - K_i'j| - 1 = -K_ij K_i'j on I.
Exact integer arithmetic for (a)-(c); floats only for the statistic in (d)."""
import sys, math, random, itertools
import numpy as np
from math import comb

def chi_table(q):
    sq = set((x * x) % q for x in range(1, q))
    return [0] + [1 if x in sq else -1 for x in range(1, q)]

def paley(q):
    ch = chi_table(q)
    return [[ch[(j - i) % q] for j in range(q)] for i in range(q)]

def check_ab(q, rho=0.995, seed=1):
    K = paley(q)
    A = [i for i in range(q) if K[i][0] == 1]
    assert len(A) == (q - 1) // 2 and 0 not in A
    # column sign counts on A, columns v != 0
    minpos = minneg = q
    for v in range(1, q):
        pos = sum(1 for i in A if K[i][v] == 1)
        neg = sum(1 for i in A if K[i][v] == -1)
        zer = sum(1 for i in A if K[i][v] == 0)
        assert zer == (1 if v in A else 0)
        assert pos - neg == (-1 - K[0][v]) // 2   # derived: sum_A K_iv = (-1-chi(v))/2
        minpos = min(minpos, pos); minneg = min(minneg, neg)
    # random selection of m rows from A
    rng = random.Random(seed)
    m = int(round(rho * (q - 1) / 2))
    S = sorted(rng.sample(A, m))
    deleted = len(A) - m
    ok_norm = all(any(K[i][v] == 1 for i in S) and any(K[i][v] == -1 for i in S) for v in range(1, q))
    # parity of exponent 1+K_iv is singleton indicator of owner i = v
    par_ok = all(((1 + K[i][v]) % 2) == (1 if i == v else 0) for i in S for v in range(1, q))
    # threshold bits [a>=1],[a>=2]; parity = sum mod 2
    def sig(i):
        bits = []
        for v in range(1, q):
            a = 1 + K[i][v]
            bits.append(1 if a >= 1 else 0); bits.append(1 if a >= 2 else 0)
        return tuple(bits)
    sigs = {i: np.array(sig(i), dtype=np.uint8) for i in S}
    thr_par_ok = all(((sigs[i][2*(v-1)] + sigs[i][2*(v-1)+1]) % 2) == ((1 + K[i][v]) % 2) for i in S for v in range(1, q))
    sums = set()
    npairs = 0
    for i, j in itertools.combinations(S, 2):
        sums.add(((sigs[i] ^ sigs[j])).tobytes()); npairs += 1
    sidon = (len(sums) == npairs)
    # sums distinct from zero sum too (zero = i==j)
    print(f"q={q}: |A|={len(A)}, min #(+1) per column on A={minpos}, min #(-1)={minneg}, (q-1)/4={(q-1)/4:.2f};"
          f" m={m}, deleted rows={deleted}; both exponents 0,2 in every column after selection: {ok_norm};"
          f" parity=owner indicator: {par_ok}; parity=sum of threshold bits: {thr_par_ok};"
          f" Sidon (distinct unary signature sums over {npairs} pairs): {sidon}")
    return minpos, minneg

def rm_dim(q, d, c=1):
    r = math.ceil(math.log2(q)) + c   # n = 2^r labels, n = O(s)
    dim = sum(comb(r, j) for j in range(d + 1))
    return r, dim

def margin_stat(q, trials=3, seed=2, kind="uniform"):
    K = np.array(paley(q), dtype=np.float64)
    A = [i for i in range(q) if K[i][0] == 1]
    KA = K[np.ix_(A, list(range(1, q)))]          # rows A, columns I
    n = q - 1
    rng = np.random.default_rng(seed)
    out = []
    for t in range(trials):
        if kind == "uniform":
            x = rng.random(n)                      # values in [0,1]; sigma = 1 after scaling
        else:
            x = np.zeros(n); x[: n // 2] = 1.0     # extreme two-point profile
        xbar = x.mean()
        perm = rng.permutation(n)
        xp = x[perm]
        G = -(KA * xp) @ KA.T                      # G[i,i'] = sum_j c_j x_{pi(j)}
        iu = np.triu_indices(len(A), 1)
        dev = G[iu] - 2 * xbar
        out.append((dev.min(), dev.max(), np.abs(dev).max()))
    scale = math.sqrt(n * math.log(n))
    for (mn, mx, ab) in out:
        print(f"  q={q} [{kind}] sigma=1: min dev={mn:.3f}, max dev={mx:.3f}, max|dev|={ab:.3f},"
              f" sqrt(n log n)={scale:.1f}, ratio={ab/scale:.3f};"
              f" with sigma=q^-3/4: max|dev|={ab*q**-0.75:.4f}")

if __name__ == "__main__":
    for q in [43, 67, 131, 227, 383]:
        check_ab(q)
    for d in [1, 2, 3, 4]:
        for q in [10**3, 10**4, 10**6, 10**9]:
            r, dim = rm_dim(q, d)
            print(f"RM: d={d}, s~{q}, r={r}, dim RM(d,r)={dim}, m~s/2={q//2}, dim<m: {dim < q//2}")
    for q in [227, 1019, 2003, 4003]:
        margin_stat(q, kind="uniform")
        margin_stat(q, trials=1, kind="twopoint")
