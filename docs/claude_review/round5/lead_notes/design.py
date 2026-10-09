# (1) Residue design x -> x mod p_q: max weighted pair demand over moduli q^a with (q-1)q^(a-1) < 2M.
# (2) Multi-block Sylvester with random row permutations: exhaustive min of sum_i E_i over
#     all +-1 characters of support 4 (and support <= 4 overall).
import itertools, math, random, sys
import numpy as np
def primes_upto(n):
    s = bytearray([1])*(n+1); s[0]=s[1]=0
    for i in range(2,int(n**0.5)+1):
        if s[i]: s[i*i::i]=bytearray(len(s[i*i::i]))
    return [i for i in range(n+1) if s[i]]
def design_demand(M):
    P = primes_upto(4*M)
    odd = [q for q in P if q > 2]
    # group order m_q: inert q=3 mod 4 -> q+1 ; split q=1 mod 4 (not dividing N) -> q-1
    used = set(); assign = {}
    for q in odd:
        m = q + 1 if q % 4 == 3 else q - 1
        if m >= 2*M: continue   # pigeonhole-free moduli: make classes injective, no collisions
        cands = [p for p in P if p <= m and p not in used]
        if not cands: raise SystemExit(f"no prime for q={q}")
        p = max(cands); used.add(p); assign[q] = (p, m)
    worst = 0; arg = None
    for d in range(1, M):
        dem = 0.0
        for q, (p, m) in assign.items():
            a = 1
            mod = p
            while mod <= M and d % mod == 0 and (q - 1) * q**(a-1) < 2*M:
                dem += 2*math.log(q); a += 1; mod = p * q**(a-1)
        if dem > worst: worst, arg = dem, d
    return worst, arg, len(assign)
def sylvester(m):
    H = np.array([[1]])
    for _ in range(m): H = np.block([[H, H], [H, -H]])
    return H
def multiblock_min_support4(m, t, seed):
    H = sylvester(m); M = H.shape[0]; Hp = H[:,1:]
    rng = np.random.default_rng(seed)
    perms = [rng.permutation(M) for _ in range(t)]
    blocks = [Hp[p] for p in perms]          # row x of block i is Sylvester row perms[i][x]
    best = None; cnt_small = 0
    quads = list(itertools.combinations(range(M), 4))
    Q = np.array(quads)
    for signs in [(1,1,-1,-1),(1,-1,1,-1),(1,-1,-1,1)]:
        s = np.array(signs)
        tot = np.zeros(len(Q), dtype=np.int64)
        for B in blocks:
            T = (B[Q] * s[None,:,None]).sum(axis=1)        # (#quads, M-1)
            tot += np.abs(T).sum(axis=1) - (M-1)
        mn = tot.min(); cnt_small += int((tot < t*M//4).sum())
        best = mn if best is None else min(best, mn)
    return M, t, int(best), cnt_small
if __name__ == "__main__":
    for M in (64, 256, 1024, 4096):
        w, d, k = design_demand(M)
        print(f"design M={M}: moduli used {k}, max pair demand {w:.2f} at |x-y|={d}; 2 log M = {2*math.log(M):.2f}; ratio {w/math.log(M):.2f} log M")
    for (m, t) in ((5, 1), (5, 3), (6, 1), (6, 3), (6, 5)):
        for seed in (1, 2):
            print("multiblock", multiblock_min_support4(m, t, seed), "(M, t, min sum E over support-4 chars, #chars with sum E < tM/4)")
