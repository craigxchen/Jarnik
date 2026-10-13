# Referee's own random sample at n = 5 (d = 4): random full-rank sign profiles (own seed), T <= Tmax, L <= 3;
# admissible = (i),(ii'),(iii),(v); exact rho; zonotope witness search for sum L_i = rho.
import sys, itertools, random, time
import flint
from rzono import *
n = 5; base = 4
seed = int(sys.argv[1]); nprof = int(sys.argv[2]); Tmax = int(sys.argv[3])
random.seed(seed)
pats = []
for r_ in range(1, n):
    for S in itertools.combinations(range(n), r_):
        S = frozenset(S); Sc = frozenset(range(n)) - S
        if min(tuple(sorted(S)), tuple(sorted(Sc))) == tuple(sorted(S)): pats.append(S)
def full_rank(prof):
    others = [x for x in range(n) if x != base]
    rows = [[(1 if x in S else 0) - (1 if base in S else 0) for x in others] for S in prof]
    return flint.fmpz_mat(rows).rank() == 4
profiles = []
while len(profiles) < nprof:
    r = random.randint(n, 3 * n)
    prof = list({random.choice(pats) for _ in range(r)})
    if full_rank(prof): profiles.append(prof)
t0 = time.time()
seen = {}
for prof in profiles:
    for T in range(1, Tmax + 1):
        for L in range(1, 4):
            others, cons = functionals(n, prof, T, L, base)
            A = points(4, cons, T)
            if len(A) == (2 * T + 1) ** 4: continue
            key = tuple(A)
            if key in seen: continue
            rh = rho(A, T, exact=False)
            rngs = []
            for S in prof:
                lam = [(1 if x in S else 0) - (1 if base in S else 0) for x in others]
                vals = [sum(a * b for a, b in zip(lam, k)) for k in A]
                rngs.append(max(vals) - min(vals))
            if max(rngs) + 1 > rh:
                seen[key] = None; continue
            rh = rho(A, T, exact=True)
            w = find_zonotope(cons, 4, rh, E=2, node_limit=int(sys.argv[4]) if len(sys.argv) > 4 else 3_000_000)
            seen[key] = (rh, w is not None)
            print(f"profile {[sorted(S) for S in prof]} T={T} L={L}: |A|={len(A)} rho={rh} max range={max(rngs)} "
                  f"witness zeta=rho: {'FOUND' if w else 'not found (search limit)'} ({time.time()-t0:.0f}s)", flush=True)
adm = [v for v in seen.values() if v]
print(f"distinct admissible A: {len(adm)}; witness found for {sum(1 for v in adm if v[1])}")
