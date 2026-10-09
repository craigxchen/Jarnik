# brute-force search for rational points of the C5-ansatz variety: cubic in rho12
from c5ansatz import *
import random, sys, time
random.seed(int(sys.argv[1]) if len(sys.argv) > 1 else 0)
R = int(sys.argv[2]) if len(sys.argv) > 2 else 6
NT = int(sys.argv[3]) if len(sys.argv) > 3 else 2000

def interp(xs, ys):
    n = len(xs); res = []
    for i in range(n):
        num = [Fr(1)]; den = Fr(1)
        for j in range(n):
            if j != i:
                num = mul(num, linpoly(1, -xs[j])); den *= xs[i] - xs[j]
        res = add(res, scal(ys[i] / den, num))
    return res

hits = []
t0 = time.time()
for trial in range(NT):
    vals = random.sample(range(-R, R + 1), 9)
    rho = {e: Fr(v) for e, v in zip(EDGES[1:], vals[:4])}
    sigma = {T: Fr(v) for T, v in zip(TRIPLES, vals[4:])}
    bad = set(rho.values()) | set(sigma.values())
    def N(x):
        r = dict(rho); r[(1, 2)] = x
        D = prod([linpoly(-1, sigma[T]) for T in triples_of(1)] + [linpoly(-1, sigma[T]) for T in triples_of(2)])
        return det_frac(matrix(r, sigma)) * ev(D, x)
    xs = [Fr(10**6 + 17 * k, 1) for k in range(5)]
    p = interp(xs, [N(x) for x in xs])
    if not p:
        continue
    q, rem = divmod_(p, linpoly(1, -sigma[(1, 2, 4)]))
    assert not rem
    for r in rational_roots(q):
        if r in bad:
            continue
        rr = dict(rho); rr[(1, 2)] = r
        sol = solve_curve(rr, sigma)
        if sol is None:
            continue
        hits.append((rr, sigma, sol))
        print('HIT', {e: str(v) for e, v in rr.items()}, {T: str(v) for T, v in sigma.items()}, flush=True)
print('trials', NT, 'hits', len(hits), 'time', round(time.time() - t0, 1))
