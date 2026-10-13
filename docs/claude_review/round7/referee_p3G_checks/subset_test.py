"""Does Theorem 4.1 transfer to members using a SUBSET of a verified prime set?
p3G.md 4.2 says Thm 4.1 'covers every b-copy Paley profile with b(M-1) <= r built on these primes'.
Such a member has the subset as its prime set; its closure uses only relations supported on the subset.
usage: python3 subset_test.py M b K seed rsub nsamples
"""
import sys, json, random
from rcore import rand_primes_like_author, gauss, Lprime, tight_list
from rclosure import closure, status, check

M, b, K, seed, rsub, ns = int(sys.argv[1]), int(sys.argv[2]), float(sys.argv[3]), int(sys.argv[4]), int(sys.argv[5]), int(sys.argv[6])
Kv = int(K) if K == int(K) else K
P = int(M ** Kv)
ps = rand_primes_like_author(b * (M - 1), P, 2 * P, random.Random(seed))
X = 2 * M
L = Lprime(X, set(ps))
srng = random.Random(12345)
for t in range(ns):
    sub = sorted(ps[:rsub]) if t == 0 else sorted(srng.sample(ps, rsub))
    gp = [gauss(p) for p in sub]
    rels = tight_list(gp, L)
    Q, log = closure(len(gp), rels)
    ok, exact, half = check(gp, L, rels, log)
    st = status(Q, len(gp))
    print(json.dumps(dict(base=f'M={M} b={b} K={sys.argv[3]} seed={seed}', subset=('first' if t == 0 else 'random'),
                          rsub=rsub, tight=len(rels), certs=ok, free=st.count(0), half_pinned=len(half),
                          exact_pinned=len(exact))), flush=True)
