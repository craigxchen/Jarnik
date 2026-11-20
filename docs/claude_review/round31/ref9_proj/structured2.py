import sys; sys.path.insert(0, '/tmp/claude-0/-home-user-Jarnik/ba4929a4-8665-5790-985a-3f09130a14b9/scratchpad/round31/ref9_proj')
from structured import *
from statistics import median
P1 = [p for p in primerange(5, 400) if p % 4 == 1]
def run(N, C, MM, lab):
    W = windows(N, C, MM)
    seen = set(); out = []
    for cl in W:
        g, zs, Np = normalise(cl)
        key = min(tuple(sorted(mul(u, z) if not cj else conj(mul(u, z)) for z in zs)) for u in UNITS for cj in (0, 1))
        if key in seen: continue
        seen.add(key)
        M = len(zs); primes, pis, e, rows = rows_of(zs, Np)
        best = None; bestpair = None
        cs = [[1 if t in (i, j) else 0 for t in range(M)] for i, j in combinations(range(M), 2)]
        cs += [[2 if t == i else 0 for t in range(M)] for i in range(M)]
        if M <= 6:
            cs += [list(c) for c in product(range(-1, 3), repeat=M) if sum(c) == 2 and 3 <= sum(map(abs, c)) <= 4]
        for c in cs:
            T, b, u, v = certificate(zs, Np, primes, pis, e, rows, c)
            l1 = sum(map(abs, c))
            q = log((1+l1)**2*T)/log(Np)
            if best is None or q < best: best = q
            if l1 == 2 and max(c) == 1:
                qp = log(T)/log(Np) if T > 1 else 0.0
                bestpair = qp if bestpair is None else min(bestpair, qp)
        out.append((M, best, bestpair))
    if out:
        bs = sorted(o[1] for o in out); bp = sorted(o[2] for o in out)
        print(f"{lab} N'={N:.3e} C={C} M>={MM}: {len(out)} windows (mod symmetry), Mmax={max(o[0] for o in out)}; "
              f"log((1+|c|)^2 T)/log N': min {bs[0]:.3f} median {median(bs):.3f}; pair log T/log N': min {bp[0]:.3f} median {median(bp):.3f}; "
              f"fraction with a pair T=1: {sum(o[2]==0 for o in out)/len(out):.3f}", flush=True)
    else:
        print(f"{lab} C={C} M>={MM}: none", flush=True)
for k in (6, 7, 8, 9, 10):
    N = 1
    for p in P1[:k]: N *= p
    for C, MM in ((2, 3), (4, 4), (8, 5)):
        run(N, C, MM, f"first{k}")
# random squarefree products of 8 primes in [100,1000]
import random
random.seed(1)
P2 = [p for p in primerange(100, 1000) if p % 4 == 1]
for t in range(6):
    ps = random.sample(P2, 8); N = 1
    for p in ps: N *= p
    for C, MM in ((4, 3), (8, 4)):
        run(N, C, MM, f"rand8#{t}")
