# comp(F_A, mu) on random integer Mobius curves passing the transversality filter (class mu verified).
import sys, random; sys.path.insert(0, '.')
from reflib import *
cube = {(i,j,k): meet3(i,j,k) for i in (1,2) for j in (3,4) for k in (5,6)}
Ekeys = [(1,3,5),(1,3,6),(1,4,5),(1,4,6),(2,3,5),(2,4,5)]
webm = web([cube[k] for k in Ekeys])
seed = int(sys.argv[1]) if len(sys.argv) > 1 else 1
rng = random.Random(seed)
res = []
tries = 0
while len(res) < int(sys.argv[2] if len(sys.argv) > 2 else 4):
    tries += 1
    V = [[P([rng.randint(-9, 9), rng.randint(-9, 9)]), P([rng.randint(-9, 9), rng.randint(-9, 9)])] for a in range(7)]
    cls, good, info = class_and_goodness(V)
    if not good:
        continue
    if not (all(cls.get(frozenset(k), 0) == 2 for k in combinations(range(7), 2)) and sum(cls.values()) == 42):
        continue
    Br, dBr, _ = FA_brackets(V, webm)
    c, nk = comp_count(Br, dBr)
    res.append(c)
    print('mu instance', len(res), 'comp_A =', c, ' nonzero W:', nk, flush=True)
print('seed', seed, 'values', res, 'tries', tries)
