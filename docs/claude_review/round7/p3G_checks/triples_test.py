import random, time, math, json, sys
from collections import Counter
from gtight import *
from verify_certs import verify
from trend import rand_primes
out = open('triples_out.jsonl', 'a')
M = int(sys.argv[1]); b = 1
for K in [float(x) for x in sys.argv[2].split(',')]:
    for seed in range(int(sys.argv[3])):
        rng = random.Random(seed)
        X = 2 * M; r = b * (M - 1); P = int(M ** K)
        ps = rand_primes(r, P, 2 * P, rng)
        gp = []
        for p in ps:
            a, bb = two_squares(p)
            if rng.random() < 0.5: bb = -bb
            gp.append((a, bb))
        t = time.time()
        rels = tight_relations(gp, X)
        P2, log2 = closure(r, rels)
        f2 = Counter(pinned_status(P2, r))[0]
        rels3 = rels + triple_relations(gp, X)
        P3, log3 = closure(r, rels3)
        f3 = Counter(pinned_status(P3, r))[0]
        ex, hf = verify(gp, X, rels3, log3)
        nq = sum(1 for x in rels3 if x[0] == 'Q')
        d = dict(M=M, K=K, seed=seed, r=r, pair_deg=2 * sum(1 for x in rels if x[0] == 'P') / r, triple_per_col=3 * nq / r,
                 free_pairs_only=f2, free_with_triples=f3, verified_half=hf, secs=round(time.time() - t, 1))
        out.write(json.dumps(d) + '\n'); out.flush()
        print(d, flush=True)
