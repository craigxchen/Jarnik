"""Lattice (CRT) census of balanced anchored configurations with prescribed blocks (Data + heuristic).

For k = m+1 and a scale A = e^w: assign to each visible T a distinct prime n_T in [A, 2A] (random
bijection).  The anchored configurations P_0 = (1,0), P_i = (X_i, 1), X_1 = 0, realising these blocks
(n_T | X_i - X_j for all i, j in T) form the lattice
    Lambda(n) = {X in Z^{m-1} : X_i = X_j mod n_T  (i, j in T)},   index prod_T n_T^{|T|-1}.
For X in Lambda(n) with distinct coordinates, the datum with these n_T has Y_j = 1 and
    t_ij = (X_i - X_j)/G_ij,  G_ij = prod_{T visible, T contains i,j} n_T,
so H_Sigma(X) = sum_{i<j} log|X_i - X_j| - sum_{i<j} log G_ij  (exact integers; logs for output).
The blocks are pairwise coprime primes, so the two-chain condition holds, and the balance is
eta <= log 2 / w.  We LLL-reduce a basis of Lambda(n) and evaluate H on all combinations
sum c_l b_l with |c_l| <= R; the minimum is an UPPER bound for min_{X in Lambda(n)} H_Sigma(X), hence for
the census minimum at scale w.  Gaussian heuristic for comparison:
    lambda_1 ~ (index)^{1/(m-1)},   H_typ ~ C(m,2) * (log index/(m-1) - b w),  b = 2^{m-2} - 1.
Usage: python3 lattice_census.py k A samples seed
"""
import sys, random, time, json
from itertools import product, combinations
from math import log, comb
from flint import fmpz_mat
from sympy import primerange


def visible(m):
    return [frozenset(T) for s in range(2, m) for T in combinations(range(1, m + 1), s)]


def lattice_basis(m, blocks):
    """blocks: dict T -> prime.  Coordinates X_2..X_m (X_1 = 0).  Returns list of basis vectors."""
    d = m - 1
    B = [[1 if r == c else 0 for r in range(d)] for c in range(d)]  # columns as lists (vectors)

    def coord(vec, i):  # value of X_i in a vector of coordinates (X_2..X_m)
        return 0 if i == 1 else vec[i - 2]

    for T, p in blocks.items():
        idx = sorted(T)
        i0 = idx[0]
        for i in idx[1:]:
            # condition X_i - X_{i0} = 0 mod p
            v = [(coord(b, i) - coord(b, i0)) % p for b in B]
            js = [j for j in range(d) if v[j]]
            if not js:
                continue
            j = js[0]
            inv = pow(v[j], -1, p)
            for l in range(d):
                if l != j and v[l]:
                    f = (v[l] * inv) % p
                    B[l] = [B[l][r] - f * B[j][r] for r in range(d)]
            B[j] = [p * x for x in B[j]]
        # occasional reduction keeps entries small
    M = fmpz_mat(B)
    M = M.lll()
    return [[int(M[r, c]) for c in range(d)] for r in range(d)]


def H_of(X, m, G):
    full = [0] + X
    if len(set(full)) < m:
        return None
    h = 0.0
    for i in range(m):
        for j in range(i + 1, m):
            h += log(abs(full[i] - full[j])) - G[(i + 1, j + 1)]
    return h


def run(k, A, samples, seed, R=2):
    m = k - 1
    rng = random.Random(seed)
    Ts = visible(m)
    V = len(Ts)
    primes = list(primerange(A, 2 * A))
    if len(primes) < V:
        return {"k": k, "A": A, "error": "not enough primes in [A,2A]"}
    w = log(A)
    b = 2 ** (m - 2) - 1
    res = []
    t0 = time.time()
    for s in range(samples):
        chosen = rng.sample(primes, V)
        blocks = dict(zip(Ts, chosen))
        G = {}
        for i in range(1, m + 1):
            for j in range(i + 1, m + 1):
                G[(i, j)] = sum(log(p) for T, p in blocks.items() if i in T and j in T)
        basis = lattice_basis(m, blocks)
        best = None
        for c in product(range(-R, R + 1), repeat=m - 1):
            if not any(c):
                continue
            X = [sum(c[l] * basis[l][r] for l in range(m - 1)) for r in range(m - 1)]
            h = H_of(X, m, G)
            if h is not None and (best is None or h < best):
                best = h
        logindex = sum((len(T) - 1) * log(p) for T, p in blocks.items())
        res.append((best, logindex))
    hs = [r[0] for r in res if r[0] is not None]
    logindex = sum(r[1] for r in res) / len(res)
    heur = comb(m, 2) * (logindex / (m - 1) - b * w)
    return {"k": k, "A": A, "w": round(w, 3), "samples": samples, "minH": round(min(hs), 2),
            "meanH": round(sum(hs) / len(hs), 2), "minH/w": round(min(hs) / w, 2),
            "meanH/w": round(sum(hs) / len(hs) / w, 2), "heuristic_H_typ/w": round(heur / w, 2),
            "sec": round(time.time() - t0, 1)}


if __name__ == "__main__":
    k, A, samples, seed = map(int, sys.argv[1:5])
    print(json.dumps(run(k, A, samples, seed)), flush=True)
