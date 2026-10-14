import itertools, random, sys
from tclib import *
random.seed(1)
def formula(u, s, L):
    supp = sum(1 for x in u if x != 0)
    from collections import Counter
    t = s - max(Counter(u).values())
    return min(2*supp, 2*L + 2*t, s + L)
bad = 0; tot = 0
for s in range(1, 7):
    mons = cube_monomials(s)
    for L in range(0, s):
        pts = box_points(s, [(list(range(s)), L)])
        # all value patterns: u in {0,1,2,3,...}^s up to s values; test patterns via set partitions with small values
        tested = set()
        for u in itertools.product(range(-1, 3), repeat=s) if s <= 5 else [tuple(random.choice([-1,0,1,2,3,5]) for _ in range(s)) for _ in range(300)]:
            if all(x == 0 for x in u): continue
            nu, rho = nu_mod(pts, u, mons=mons)
            f = formula(u, s, L)
            tot += 1
            if nu != f:
                bad += 1
                if bad < 30: print('MISMATCH s=%d L=%d u=%s nu=%d formula=%d rho=%d' % (s, L, u, nu, f, rho))
        print('s=%d L=%d done rho=%d' % (s, L, rho), flush=True)
print('total', tot, 'bad', bad)
