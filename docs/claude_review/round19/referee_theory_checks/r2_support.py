"""Referee check r2 (independent). N(u) = #{T subset [7], 2<=|T|<=6 : e(T) != 0}, e(T) = sum_{i<j in T} e_ij,
for units u = prod (x_i-x_j)^{e_ij}, sum e_ij = 0, e != 0.
 (a) exhaustive over all e in {-1,0,1}^21 with sum 0 and support <= 6;
 (b) exhaustive over e in {-2..2}^21, sum 0, support <= 3;
 (c) 300000 random e (entries in -3..3, random supports), sum forced 0;
 (d) min N split by k(u) = |I| + [some vertex sum != 0] (I = vertices touched by supp e);
 (e) cross-ratio N.
Also checks q(1_T) = 0 for |T|<=1 and T = [7] (so all nonzeros are among the 119 divisors).
"""
import numpy as np, itertools, random, sys
pairs = list(itertools.combinations(range(7), 2))
subs = [T for r in range(0, 8) for T in itertools.combinations(range(7), r)]
M = np.array([[1 if (i in T and j in T) else 0 for (i, j) in pairs] for T in subs], dtype=np.int64)
valid = np.array([2 <= len(T) <= 6 for T in subs])

def N_of(E):  # E: (m,21)
    V = E @ M.T
    assert not np.any(V[:, ~valid]), "nonzero off the 119"
    return (V != 0).sum(axis=1)

def kval(e):
    I = set(); deg = [0] * 7
    for idx, (i, j) in enumerate(pairs):
        if e[idx]:
            I |= {i, j}; deg[i] += e[idx]; deg[j] += e[idx]
    return len(I) + (1 if any(deg) else 0)

stats = {}
def record(E):
    Ns = N_of(E)
    for e, n in zip(E, Ns):
        k = kval(e)
        if k not in stats or n < stats[k][0]:
            stats[k] = (int(n), tuple(int(x) for x in e))
    return int(Ns.min())

overall = 10**9
# (a)
cnt = 0
for s in range(2, 7):
    batch = []
    for supp in itertools.combinations(range(21), s):
        for signs in itertools.product((-1, 1), repeat=s):
            if sum(signs) != 0:
                continue
            e = [0] * 21
            for p, sg in zip(supp, signs):
                e[p] = sg
            batch.append(e)
            if len(batch) >= 200000:
                overall = min(overall, record(np.array(batch, dtype=np.int64))); cnt += len(batch); batch = []
    if batch:
        overall = min(overall, record(np.array(batch, dtype=np.int64))); cnt += len(batch)
print("(a) {-1,0,1}, support<=6, sum 0: vectors", cnt, " min N", overall, flush=True)
# (b)
cnt2 = 0; mn2 = 10**9; batch = []
for s in range(2, 4):
    for supp in itertools.combinations(range(21), s):
        for vals in itertools.product((-2, -1, 1, 2), repeat=s):
            if sum(vals) != 0:
                continue
            e = [0] * 21
            for p, v in zip(supp, vals):
                e[p] = v
            batch.append(e)
mn2 = record(np.array(batch, dtype=np.int64)); cnt2 = len(batch)
print("(b) {-2..2}, support<=3, sum 0: vectors", cnt2, " min N", mn2, flush=True)
# (c)
rng = random.Random(20261028)
batch = []
for _ in range(300000):
    s = rng.randint(2, 21)
    supp = rng.sample(range(21), s)
    e = [0] * 21
    for p in supp:
        e[p] = rng.choice((-3, -2, -1, 1, 2, 3))
    tot = sum(e)
    e[supp[0]] -= tot
    if not any(e):
        continue
    batch.append(e)
mn3 = record(np.array(batch, dtype=np.int64))
print("(c) random vectors", len(batch), " min N", mn3, flush=True)
print("(d) evidence: smallest N found per k(u):", {k: v[0] for k, v in sorted(stats.items())})
for k, v in sorted(stats.items()):
    print("     k =", k, "example e:", {pairs[i]: x for i, x in enumerate(v[1]) if x})
# (e)
e = [0] * 21; e[pairs.index((0, 1))] = 1; e[pairs.index((0, 2))] = -1
print("(e) cross-ratio (x1-x2)/(x1-x3): N =", int(N_of(np.array([e]))[0]), " k =", kval(e))
