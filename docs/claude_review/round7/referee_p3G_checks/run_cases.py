"""Independent re-run of p3G.md Theorem 4.1 (and partial cases) with referee code.
usage: python3 run_cases.py M b K seed [X]
Prime sets: same generator description (random.Random(seed), randrange(P,2P)|1, p=1 mod 4),
so the prime set should coincide with the author's; pmin/pmax are printed for comparison.
"""
import sys, json, time, random
from collections import Counter
from rcore import rand_primes_like_author, gauss, Lprime, tight_list
from rclosure import closure, status, check

def run(ps, X, name):
    t = time.time()
    gp = [gauss(p) for p in ps]
    r = len(gp)
    L = Lprime(X, set(ps))
    rels = tight_list(gp, L)
    kinds = Counter(('S' if len(v) == 1 and v[0][1] == 1 else 'T' if len(v) == 1 else 'P') for v, *_ in rels)
    Q, log = closure(r, rels)
    st = Counter(status(Q, r))
    ok, exact, half = check(gp, L, rels, log)
    d = dict(name=name, r=r, X=X, pmin=min(ps), pmax=max(ps), singles=kinds['S'], doubles=kinds['T'],
             pairs=kinds['P'], pair_degree=round(2 * kinds['P'] / r, 2), certs=len(log), certs_checked=ok,
             uf_free=st[0], uf_half=st[1], uf_exact=st[2], chk_exact=len(exact), chk_half=len(half),
             collapse=(len(half) == r), secs=round(time.time() - t, 1))
    print(json.dumps(d), flush=True)
    return d

if __name__ == '__main__':
    M, b, K, seed = int(sys.argv[1]), int(sys.argv[2]), float(sys.argv[3]), int(sys.argv[4])
    X = int(sys.argv[5]) if len(sys.argv) > 5 else 2 * M
    Kv = int(K) if K == int(K) else K
    P = int(M ** Kv)
    rng = random.Random(seed)
    ps = rand_primes_like_author(b * (M - 1), P, 2 * P, rng)
    run(ps, X, f'random M={M} b={b} K={sys.argv[3]} seed={seed}')
