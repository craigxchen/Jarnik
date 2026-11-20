"""Check: for the Legendre profile a_ij = 1+chi(j-i), e_j = 2 (j != 0), every integer c with sum 2 on rows A
has K^T c nonzero on at least (q-1)/2 of the columns j != 0. Referee's own code."""
import random, itertools
def chi(x, q):
    x %= q
    if x == 0: return 0
    return 1 if pow(x, (q-1)//2, q) == 1 else -1
for q in (7, 11, 19, 23, 31, 43):
    A = [i for i in range(q) if chi(0-i, q) == 1]   # rows with +1 on column 0
    worst = q
    trials = 0
    rng = random.Random(q)
    cand = []
    for i, j in itertools.combinations(A, 2): cand.append({i: 1, j: 1})
    for i in A: cand.append({i: 2})
    for _ in range(20000):
        k = rng.randint(1, len(A)); rows = rng.sample(A, k)
        c = {r: rng.randint(-3, 3) for r in rows}
        s = sum(c.values()); r0 = rng.choice(rows); c[r0] += 2 - s
        cand.append(c)
    for c in cand:
        assert sum(c.values()) == 2
        nz = sum(1 for j in range(1, q) if sum(ci*chi(j-i, q) for i, ci in c.items()) != 0)
        worst = min(worst, nz)
    print(f"q={q}: |A|={len(A)}, tested {len(cand)} c, min #nonzero v_j on j!=0: {worst}, (q-1)/2={(q-1)//2}")
