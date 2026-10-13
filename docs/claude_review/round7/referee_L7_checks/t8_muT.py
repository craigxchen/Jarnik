# Evidence checks: comp_A(mu_T) - comp_A(mu) for all 35 triples T (L7 Sec 5.4 'Equality' item 1), exact, 2 instances each;
# and lambda_{f_i}(tau gamma_A) value set (Cor 4.5 bullet), via Lemma 4.1 formula.
import sys, random; sys.path.insert(0, '.')
from reflib import *
cube = {(i,j,k): meet3(i,j,k) for i in (1,2) for j in (3,4) for k in (5,6)}
Ekeys = [(1,3,5),(1,3,6),(1,4,5),(1,4,6),(2,3,5),(2,4,5)]
webm = web([cube[k] for k in Ekeys])
rng = random.Random(77)
def rand_mob():
    return [[P([rng.randint(-9, 9), rng.randint(-9, 9)]), P([rng.randint(-9, 9), rng.randint(-9, 9)])] for a in range(7)]
def comp_of(V):
    Br, dBr, _ = FA_brackets(V, webm); return comp_count(Br, dBr)[0]
def is_mu(cls):
    return all(cls.get(frozenset(k), 0) == 2 for k in combinations(range(7), 2)) and sum(cls.values()) == 42
mus = []
while len(mus) < 2:
    V = rand_mob(); cls, good, _ = class_and_goodness(V)
    if good and is_mu(cls): mus.append(comp_of(V))
print('comp(mu):', mus)
delta = {}
for T in combinations(range(7), 3):
    vals = []
    while len(vals) < 2:
        V = rand_mob()
        # force a triple event of T at s=0: set constant terms of labels in T equal to that of T[0] (times random scalar)
        a0 = T[0]
        for a in T[1:]:
            lam = rng.choice([1, -1, 2])
            V[a] = [P([int(V[a0][0].coeffs()[0] if V[a0][0] != 0 else 0) * lam, rng.randint(-9, 9)]),
                    P([int(V[a0][1].coeffs()[0] if V[a0][1] != 0 else 0) * lam, rng.randint(-9, 9)])]
        cls, good, _ = class_and_goodness(V)
        exp = {frozenset(k): (1 if set(k) <= set(T) else 2) for k in combinations(range(7), 2)}
        exp[frozenset(T)] = 1
        if not good or {k: v for k, v in cls.items()} != exp:
            continue
        vals.append(comp_of(V))
    delta[T] = [v - mus[0] for v in vals]
    print(''.join(LAB[i] + ' ' for i in T), delta[T], flush=True)
s = sum(d[0] for d in delta.values())
print('sum of deltas', s, ' => 3*comp(mu)+sum =', 3 * mus[0] + s)
