"""Minimal defect of pencils: s = sum_cuts width(S) - 8t with dim A_t(S) >= 2 (k=5 slice, coords 1..4).
Search: zonotope sumsets of 1..3 root generators, plus 1-2 extra lattice points near them, t = ceil(MW) - j."""
import random, math
from rigidity import dimA, MW, reg
from fractions import Fraction
rng = random.Random(3)
roots = [(a, b) for a in range(5) for b in range(a+1, 5)]
def vec(a, b):
    v = [0]*5; v[a] += 1; v[b] -= 1
    return tuple(v[:4])
def cutsum(S):
    return MW(S, 4) * 8  # = sum over 15 cuts of widths
best = {}
nbr = [tuple(v) for v in __import__('itertools').product([-1,0,1], repeat=4) if any(v)]
for trial in range(2500):
    k = rng.randint(0, 3)
    gens = [vec(*rng.choice(roots)) for _ in range(k)]
    Z = {tuple([0]*4)}
    for gv in gens:
        Z = Z | {tuple(z[i]+gv[i] for i in range(4)) for z in Z}
    Z = set(Z)
    for _ in range(rng.randint(1, 2)):
        z = rng.choice(sorted(Z)); d = rng.choice(nbr)
        Z.add(tuple(z[i]+d[i] for i in range(4)))
    S = sorted(Z)
    cs = cutsum(S)
    r = reg(S, 4)
    for t in range(max(1, r-1), r+1):
        if dimA(S, 4, t) >= 2:
            s = cs - 8*t
            if s not in best or len(S) < len(best[s][0]):
                best[s] = (S, t)
for s in sorted(best)[:8]:
    S, t = best[s]
    print('defect s =', s, 't =', t, '|S| =', len(S), 'dimA =', dimA(S, 4, t), 'S =', S)
