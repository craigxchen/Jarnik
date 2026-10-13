"""Computer-verified collapse instances (p3G.md Theorem 4.1).

For each case: build the prime set, compute all tight relations (singles, doubles, pairs) at X,
run the closure (gtight.closure), then re-verify every certificate independently (verify_certs.verify).
A case is a COLLAPSE instance iff the verifier finds 2 e_k in the lattice for every column k.
"""
import json, math, random, sys, time
from collections import Counter
from gtight import is_prime, two_squares, tight_relations, closure, pinned_status
from verify_certs import verify
from trend import rand_primes

def consecutive_primes_1mod4(start, r):
    out = []; n = start + (1 - start) % 4
    while len(out) < r:
        if is_prime(n): out.append(n)
        n += 4
    return out

def gaussians(ps, rng=None):
    gp = []
    for p in ps:
        a, b = two_squares(p)
        if rng is not None and rng.random() < 0.5: b = -b
        gp.append((a, b))
    return gp

def case(name, ps, X, rng=None):
    t = time.time()
    gp = gaussians(ps, rng)
    r = len(gp)
    rels = tight_relations(gp, X)
    Pn, log = closure(r, rels)
    st = Counter(pinned_status(Pn, r))
    ex, hf = verify(gp, X, rels, log)
    kinds = Counter(x[0] for x in rels)
    d = dict(name=name, r=r, X=X, pmin=min(ps), pmax=max(ps), tight_singles=kinds['S'], tight_doubles=kinds['T'],
             tight_pairs=kinds['P'], pair_degree=2 * kinds['P'] / r, certs=len(log), verified_exact=ex,
             verified_half=hf, collapse=(hf == r), secs=round(time.time() - t, 1))
    print(json.dumps(d), flush=True)
    return d

if __name__ == '__main__':
    out = open('collapse_cases_out.jsonl', 'w')
    # (a) prime set of the round-6 P3-char fake (sharp.md Thm 5.1, A = 1, b = 33) at M = 44
    M = 44; b = 33; r = b * (M - 1)
    Z = int(max(3 * M, 21 * r * r) * math.log(M) ** 2)
    ps = consecutive_primes_1mod4(Z, r)
    for X in (2 * M, 3 * M):
        out.write(json.dumps(case(f'sharp-Thm5.1 primes M=44 b=33 (Z={Z})', ps, X)) + '\n')
    # (b) random primes p in [M^K, 2M^K], r = b(M-1), X = 2M
    for (M, b, K, seed) in [(44, 2, 3, 0), (84, 8, 3, 0), (164, 8, 3, 0), (164, 8, 4, 0), (324, 2, 4, 0),
                            (644, 1, 4, 0), (644, 1, 4.25, 1), (1284, 1, 4.5, 0)]:
        rng = random.Random(seed)
        ps = rand_primes(b * (M - 1), int(M ** K), 2 * int(M ** K), rng)
        out.write(json.dumps(case(f'random M={M} b={b} K={K} seed={seed}', ps, 2 * M, rng)) + '\n')
    # (c) window primes as in p3data (construct profile, r = 9M - 8), random start in [1e7, 1e9]
    for (M, seed) in [(44, 0), (84, 0)]:
        rng = random.Random(1000 + seed)
        ps = consecutive_primes_1mod4(rng.randrange(10**7, 10**9), 9 * M - 8)
        out.write(json.dumps(case(f'p3data-window M={M} seed={seed}', ps, 2 * M, rng)) + '\n')
