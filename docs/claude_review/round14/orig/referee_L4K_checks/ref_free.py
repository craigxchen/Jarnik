"""Independent exact forced/free test: q is free for norm K iff some parity-respecting c in (Z/n)^r
(c_j = nu_j mod 2 if n even) has <v,c> != 0 mod n for all 1<=||v||_1<=K.
Exhaustive DFS; c_1 reduced to divisors of n (unit orbits). Prints largest forced prime and all forced primes."""
import sys, math
import numpy as np
from rlib import *

def shift(arr, t, n):
    # returns boolean array B with B[(x + t) % n] = arr[x]
    return np.roll(arr, t % n)

def free_vector(n, nuv, K):
    r = len(nuv)
    even = (n % 2 == 0)
    # P[k] = boolean array of reachable <w,c> for ||w||_1 <= k over chosen prefix
    P0 = [np.zeros(n, bool) for _ in range(K + 1)]
    for k in range(K + 1): P0[k][0] = True
    divs = [d for d in range(1, n) if n % d == 0]
    def cands(j, P):
        # valid c_j: for t=1..K, t*c not in P[K-t] (P symmetric)
        out = []
        rng = range(n) if j > 0 else divs
        for c in rng:
            if even and (c % 2) != nuv[j]: continue
            ok = True
            for t in range(1, K + 1):
                if P[K - t][(t * c) % n]: ok = False; break
            if ok: out.append(c)
        return out
    def extend(P, c):
        Q = []
        for k in range(K + 1):
            acc = P[k].copy()
            for t in range(1, k + 1):
                acc |= shift(P[k - t], t * c, n)
                acc |= shift(P[k - t], -t * c, n)
            Q.append(acc)
        return Q
    sol = []
    def dfs(j, P):
        if j == r: return True
        for c in cands(j, P):
            sol.append(c)
            if dfs(j + 1, extend(P, c)): return True
            sol.pop()
        return False
    if dfs(0, P0): return list(sol)
    return None

def verify(n, c, K, nuv):
    r = len(c)
    if n % 2 == 0 and any(c[j] % 2 != nuv[j] for j in range(r)): return False
    for v in VK(r, K):
        if sum(x * y for x, y in zip(v, c)) % n == 0: return False
    return True

def run(ps, K, qmax):
    N = 1
    for p in ps: N *= p
    forced = []
    for q in primes_upto(qmax):
        if q == 2 or N % q == 0: continue
        n = (q - chi4(q)) // 4
        nuv = nu(ps, q)
        c = free_vector(n, nuv, K)
        if c is None: forced.append(q)
        else: assert verify(n, c, K, nuv)
    return forced

if __name__ == "__main__":
    sp = split_primes_upto(200)
    K = int(sys.argv[1]); r = int(sys.argv[2]); which = sys.argv[3]
    ps = sp[:r] if which == "first" else sp[-r:]
    h = K // 2
    B = r  # h<=1
    AK = K * (2 * (K + 1) ** 2 * B + 2 * K + 1)
    qmax = 4 * AK + 5
    f = run(ps, K, qmax)
    print(f"K={K} r={r} {which} ps={ps} qmax=4A_K+5={qmax} q*={max(f) if f else None} forced={f}")
