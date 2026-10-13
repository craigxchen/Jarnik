# Other six-point event sets (L7 Sec 5.4 'Other event sets'): exact comp(mu), comp(mu_T); value 3comp(mu)+sum delta
# (equals lambda(beta) only under the equality assumption -> evidence).
import sys, random; sys.path.insert(0, '.')
from reflib import *
cube = {(i,j,k): meet3(i,j,k) for i in (1,2) for j in (3,4) for k in (5,6)}
allk = sorted(cube)
variants = {'omit a2b1c2,a2b2c2 (=E)': [(2,3,6),(2,4,6)], 'omit a1b1c1,a2b2c1': [(1,3,5),(2,4,5)],
            'omit a1b1c1,a2b2c2': [(1,3,5),(2,4,6)], 'omit a1b2c2,a2b1c1': [(1,4,6),(2,3,5)]}
for name, om in variants.items():
    webm = web([cube[k] for k in allk if k not in om])
    rng = random.Random(9)
    def rand_mob():
        return [[P([rng.randint(-9, 9), rng.randint(-9, 9)]), P([rng.randint(-9, 9), rng.randint(-9, 9)])] for a in range(7)]
    def comp_of(V):
        Br, dBr, _ = FA_brackets(V, webm); return comp_count(Br, dBr)[0]
    while True:
        V = rand_mob(); cls, good, _ = class_and_goodness(V)
        if good and all(cls.get(frozenset(k), 0) == 2 for k in combinations(range(7), 2)) and sum(cls.values()) == 42:
            cm = comp_of(V); break
    tot = 0
    for T in combinations(range(7), 3):
        while True:
            V = rand_mob(); a0 = T[0]
            for a in T[1:]:
                lam = rng.choice([1, -1, 2])
                V[a] = [P([int(V[a0][0].coeffs()[0] if V[a0][0] != 0 else 0) * lam, rng.randint(-9, 9)]),
                        P([int(V[a0][1].coeffs()[0] if V[a0][1] != 0 else 0) * lam, rng.randint(-9, 9)])]
            cls, good, _ = class_and_goodness(V)
            exp = {frozenset(k): (1 if set(k) <= set(T) else 2) for k in combinations(range(7), 2)}
            exp[frozenset(T)] = 1
            if good and cls == exp: break
        tot += comp_of(V) - cm
    print(name, ' comp(mu) =', cm, ' 3comp(mu)+sum delta_T =', 3 * cm + tot, flush=True)
