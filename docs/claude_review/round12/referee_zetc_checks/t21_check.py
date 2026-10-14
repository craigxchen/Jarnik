# Exact test of Theorem 2.1 on A_S(L) for s<=6 (s=6 subset), all L<s.
import itertools, random, sys
from rlib import *
def formula(u, s, L):
    sig = sum(1 for x in u if x != 0)
    from collections import Counter
    t = s - max(Counter(u).values())
    return min(2*sig, 2*L + 2*t, s + L)
random.seed(7)
tot = bad = 0
for s in range(1, int(sys.argv[1])+1):
    for L in range(0, s):
        if (s, L) == (1, 0): continue
        A = [k for k in cube_points(s) if abs(sum(k)) <= L]
        E = Exact(A)
        assert E.rho == s + L, (s, L, E.rho)
        us = set()
        vals = (-2,-1,0,1,2,3) if s <= 4 else (0,1,2,3,5)
        if s <= 4:
            us |= set(itertools.product(vals, repeat=s))
        for _ in range(150):
            us.add(tuple(random.choice((0,0,1,1,2,-3,7,4,5)) for _ in range(s)))
        for _ in range(50):
            us.add(tuple(random.randint(-20,20) for _ in range(s)))
        # additive relations u_x+u_y=u_z etc.
        us.add(tuple(range(1, s+1))); us.add(tuple([1,2,3,4,5,6][:s]))
        for u in us:
            a, b = E.nu(u), formula(u, s, L)
            tot += 1
            if a != b:
                bad += 1; print("MISMATCH", s, L, u, a, b)
        print("s", s, "L", L, "|A|", len(A), "rho", E.rho, "H", E.H, "tested", len(us), flush=True)
print("total", tot, "bad", bad)
